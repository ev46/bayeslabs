package com.bayeslabs.causal.model;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.TreeMap;

public class Event {
	
	private String name;
	private int id;
	
	private int persitence = 1;
	private float continuation = 1.0f;
	private final TreeMap<Integer,Float> elicited = new TreeMap<Integer,Float>();
	private final HashMap<Integer,Float> observations = new HashMap<Integer, Float>();
	private float[] marginal = null;
	

	private ArrayList<Integer> causes;
	private ArrayList<Integer> effects;
	
	public ArrayList<Integer> getCauses() {
		return causes;
	}
	public ArrayList<Integer> getEffects() {
		return effects;
	}
	public Event(int id)
	{
		this.id = id;
		elicited.put(0, 0.0f);
		causes = new ArrayList<Integer>();
		effects = new ArrayList<Integer>();
		this.name=id+"";
	}
	public void setMarginal(float[] marginal)
	{
		this.marginal = marginal;
	}
	public float[] getMarginal()
	{
		return this.marginal;
	}
	
	public String getName() {
		return name;
	}
	public void setName(String name) {
		this.name = name;
	}
	public Event(int id,String name)
	{
		this.id = id;
		elicited.put(0, 0.0f);
		causes = new ArrayList<Integer>();
		effects = new ArrayList<Integer>();
		this.name=name;
	}
	
	/**
	 * Adds a new cause to the event, uses the default P(event | cause) = 1.0f
	 * @param id
	 */
	public void addCause(Integer id)
	{
		TreeMap<Integer,Float> temp = new TreeMap<Integer,Float>(elicited);
		elicited.clear();
		causes.add(id);
		int newkey;
		float value;
		for(int key : temp.keySet())
		{
			value = temp.get(key);
			newkey = key << 1; // take the current key and shift its bits left by one to get the new index
			elicited.put(newkey, value);
		}
		elicited.put(1, 1.0f); // last item always has Power Set index of 1, def prob is 1.0
	}
//	public void addCause(Integer id, boolean )
//	{
//		TreeMap<Integer,Float> temp = new TreeMap<Integer,Float>(elicited);
//		elicited.clear();
//		causes.add(id);
//		int newkey;
//		float value;
//		for(int key : temp.keySet())
//		{
//			value = temp.get(key);
//			newkey = key << 1; // take the current key and shift its bits left by one to get the new index
//			elicited.put(newkey, value);
//		}
//		elicited.put(1, 1.0f); // last item always has Power Set index of 1, def prob is 1.0
//	}
	/**
	 * Removes a cause and its associated combinations.
	 * @param id
	 */
	public void removeCause(int id)
	{
		int ithbit = causes.size() - causes.indexOf(id); // determine the ith bit (start with 1 not 0)
		int ei = 1 << (causes.indexOf(id) - (causes.size()-1)); // find cause's index inside elicited array
		TreeMap<Integer,Float> temp = new TreeMap<Integer,Float>(elicited);
		elicited.clear();
		int MASK = 1 << (ithbit - 1);
		int newkey;
		float value;
		// now use the index of a single cause and remove any combinations with that cause
		for(int key : temp.keySet())
		{
			if((key & MASK) == MASK)
			{
				// we do nothing with it.
			}else{
				if(key > ei)
				{
					value = temp.get(key);
					// shift bits to the right by 1
					newkey = key >> 1;
					elicited.put(newkey, value);
				}else{
					elicited.put(key, temp.get(key));
				}
			}
		}
		causes.remove(new Integer(id));
		
	}
	public int getId() {
		return id;
	}
	
	/**
	 * Method adds a new or sets an existing elicited probability, note the calling
	 * client is responsible for verifying the bits being set.  We only check for 
	 * out of bounds (greater than 2^causes.size).
	 * @param index represents index in the powerset representation of combinations of causes
	 * @param probability new probability
	 */
	public void addElicittion(int index, float probability)
	{
		if(index < (1 << causes.size()))
			elicited.put(index, probability);
		
	}

	/**
	 * Removes an elicited probability using a power set index for combinations, e.g. binary 101 is index 5.
	 * This method should only be used to remove combinations of causes and not individual causes.
	 * @param index
	 */
	public void removeElicitation(int index)
	{
		elicited.remove(index);
	}
	
	public void addEffect(int effect)
	{
		effects.add(effect);
	}
	public void removeEffect(int effect)
	{
		effects.remove(new Integer(effect));
	}
	
	public String toString()
	{
		return elicited.toString();
	}

	public float getLeak() {
		
			return elicited.get(0);
		
	}

	public void setLeak(float leak) {
		this.elicited.put(0, leak);
	}

	public int getPersitence() {
		return persitence;
	}

	public void setPersitence(int persitence) {
		this.persitence = persitence;
	}

	public float getContinuation() {
		return continuation;
	}

	public void setContinuation(float continuation) {
		this.continuation = continuation;
	}
	public TreeMap<Integer, Float> getElicited() {
		return elicited;
	}
	public void addObservation(int time, float value)
	{
		this.observations.put(time, value);
	}
	public void clearObservationAt(int time)
	{
		this.observations.remove(time);
	}
	public void clearObservations()
	{
		this.observations.clear();
	}
	public HashMap<Integer,Float> getObservations()
	{
		return observations;
	}
	public float getObservation(int time)
	{
		return observations.get(time);
	}
	
	/**
	 * Used to preview the CPT when used for user interaction with the CPT.
	 * @return
	 */
	public float[] priviewCPT()
	{
		
		int active;
		float numerator;
		float denominator;
		float[] cpt = new float[1 << causes.size()];
		
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
				for(int bit=0;bit<causes.size();bit++) // here we use causes, slightly worse worst case
				{
					MASK = 1<<bit;
					if((index & MASK) == MASK) // if bit is set
					{
						nindex = index^MASK; // if set, clear the bit
						numerator *= (1 - cpt[nindex]); // lookup probability in cpt with nindex
						for(int nbit= ((bit+1) % causes.size());;nbit++) // here we find the next set bit and clear it
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
		
		return updateWithLeak(cpt);
	}
	
	private float[] updateWithLeak(float[] cpt)
	{
		if(cpt[0] > 0)
		{
			for(int i=1; i<cpt.length; i++)
			{
				cpt[i] = (1 - (1-cpt[i])*(1-cpt[0]));
			}
		}
		return cpt;
	}
	
	private int pop(int x)
	{
		x = x - ((x>>1) & 0x55555555);
		x = (x & 0x33333333) + ((x>>2) & 0X33333333);
		x = (x + ( x>>4 )) & 0X0F0F0F0F;
		x = x + (x>>8);
		x = x + (x >>16);
		return x & 0x0000003F;
	}

}
