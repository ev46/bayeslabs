package com.bayeslabs.causal.model;

import java.util.HashMap;
import java.util.TreeMap;



public class Node {
	
	private int[] causes; // holds topo_sort sample indices of its causes
	private int persistance;  // how long the effect lasts at P(T)
	private float continuation;  // P(T+1 | T)
	private int id;
	private float[] cpt;
	private final HashMap<Integer,Float> observations; // this includes upstream and downstream evidence
	private float[] marginal; // computed marginals
	
	
	
	
	

	public Node(Event event, int[] topocauses, int steps)
	{
		this.id = event.getId();
		this.marginal = new float[steps+1]; // time 0 is just static calculations
		this.causes = topocauses;
		this.persistance = event.getPersitence();
		this.continuation = event.getContinuation();
		this.observations = new HashMap<Integer,Float>(event.getObservations());
		this.cpt = new float[1 << causes.length];
		//cpt[0] = 0f;
		this.buildCPT(event.getElicited());
		
		
	}
	public float getMarginal(int time) {
		return marginal[time];
	}
	
	public float[] getMarginal()
	{
		return marginal;
	}
	
	public int getSteps()
	{
		return marginal.length;
	}
	public void setMarginal(int time, float probability) {
		this.marginal[time] = probability;
	}
	public int[] getCauses() {
		return causes;
	}
	public float getConditional(int index)
	{
		return cpt[index];
	}
	public boolean isRoot()
	{
		if(causes.length > 0)return false; else return true;
	}
	
	private void buildCPT(TreeMap<Integer,Float> elicited)
	{
		int active;
		float numerator;
		float denominator;
		
		if(elicited.containsKey(0)){ cpt[0] = elicited.get(0);}
		
		for(int index = 1; index < cpt.length; index++)
		{
			numerator = 1f;
			denominator = 1f;
			if(elicited.containsKey(index)) // this should take care of singletons & other specified probabilities
			{
				cpt[index] = elicited.get(index);
				continue;
			}
			active = this.pop(index);// number of set bits
			if(active < 3)
			{// NOR active bits
				int lastvalue = index;
				int nextbit;
				for(int i=0;i<active;i++) // using active count
				{
					//temp = index ^ (1<<bit); // clear the bit
					//if(index != temp) //check if bit was set
					nextbit = lastvalue & -lastvalue; // find next lower set bit
					numerator *= (1 - cpt[nextbit]); // look up single probability
					lastvalue = lastvalue ^ nextbit; //clear that bit in the next iteration
				}
				cpt[index] = (1 - numerator); //set the cpt value
			}else
			{// RNOR active bits
				int nindex;
				int dindex;
				int MASK;
				for(int bit=0;bit<causes.length;bit++) // here we use causes, slightly worse worst case
				{
					MASK = 1<<bit;
					if((index & MASK) == MASK) // if bit is set
					{
						nindex = index^MASK; // if set, clear the bit
						numerator *= (1 - cpt[nindex]); // lookup probability in cpt with nindex
						for(int nbit= ((bit+1) % causes.length);;nbit++) // here we find the next set bit and clear it
						{
							MASK = 1 << nbit; // create new mask
							if((nindex & MASK) == MASK) // if bit is set
							{
								dindex = nindex^(1<<nbit); // clear it
								denominator *= (1 - cpt[dindex]); // lookup
								break; // we only need the next one
							}
						}
						
					}else{
						continue; // not active bit
					}
				}
				// here we set the RNOR prob (check for den == 0 or num > den)
				if(denominator < numerator)
				{
					cpt[index] = 0f;
				}
				else if(denominator == 0)
				{
					cpt[index] = 1f;
				}
				else {cpt[index] = (1f - (numerator/denominator));}
			}
		}// end of index iteration
		this.updateWithLeak();
	}
	
	public int getPersistance() {
		return persistance;
	}

	public float getContinuation() {
		return continuation;
	}

	public int getId() {
		return id;
	}
	
	public boolean hasObservation(int time)
	{
		return observations.containsKey(time);
	}
	public float getObservationAt(int time)
	{
		return observations.get(time);
	}

	public static int countBitsSet(int i)
	{
		i = i - ((i>>1) & 0x55555555);
		i = (i & 0x33333333) + ((i>>2) & 0x33333333);
		return ((i + (i>>4) & 0xF0F0F0F) * 0x1010101) >> 24;
	}
	
	// counts number of bits set - uses ~ 20 instructions. (Hacker's Delight, p.65)
	private int pop(int x)
	{
		x = x - ((x>>1) & 0x55555555);
		x = (x & 0x33333333) + ((x>>2) & 0X33333333);
		x = (x + ( x>>4 )) & 0X0F0F0F0F;
		x = x + (x>>8);
		x = x + (x >>16);
		return x & 0x0000003F;
	}
	
	/**
	 * If the leak probability is not zero we update the CPT to reflect that
	 */
	private void updateWithLeak()
	{
		if(cpt[0] > 0)
		{
			for(int i=1; i<cpt.length; i++)
			{
				cpt[i] = (1 - (1-cpt[i])*(1-cpt[0]));
			}
		}
	}
	
	public String printCPT()
	{
		String cpts = "";
		for(int i=0; i<cpt.length;i++)
		{
			cpts += i+"="+cpt[i]+", ";
		}
		return cpts;
	}

}
