package com.bayeslabs.causal.sampling;

import java.util.Random;
import com.bayeslabs.causal.model.Node;

public class Prediction implements Sampler {

	private Node[] nodes;
	private int number;

	
	boolean[] samplestate;
	int[] samplecount; // prior counts

	
	public Prediction(Node[] nodes, int numberofsamples)
	{
		this.nodes = nodes;
		this.number = numberofsamples;
	}
	@Override
	
	/**
	 * Complexity O(TSN)
	 */
	public void run() {
		if(nodes.length < 1)return;
		Random R = new Random();
		
		for(int time = 0; time < nodes[0].getSteps(); time++)
		{
			samplecount = new int[nodes.length];
			for(int i = 0; i<number; i++)
			{
				samplestate = new boolean[nodes.length];
				
				for(int n=0; n<nodes.length;n++)
				{
					float prob = 0f;
					if(nodes[n].isRoot())
					{// has no causes, use leak, nodes are already toposorted
						
						if(nodes[n].hasObservation(time))
						{
							prob = nodes[n].getObservationAt(time);
						}else{
							
						
							if(time==0){
								prob = nodes[n].getConditional(0);
							}
							else{
								// lookup t - 1 and multiply by continuation
								prob =  1.0f - ((1.0f -  (nodes[n].getMarginal(time - 1) * nodes[n].getContinuation())) * (1.0f - nodes[n].getConditional(0)));  
								
							}
						}
							
					}else
					{ // this is caused by something else
						
						// observations
						if(nodes[n].hasObservation(time))
						{
							prob = nodes[n].getObservationAt(time);
						}else{
						// first we create a cpt index out of current sample state
							int bit = nodes[n].getCauses().length - 1; // we are starting with high order bits
							int index = 0;
							for(int cause : nodes[n].getCauses())
							{
								if(samplestate[cause])
									index|= 1<<bit;	//set the bit
								bit--;
							}
							// check status from last time
							prob = nodes[n].getConditional(index);
						}
					}
					
					
					if(prob >= R.nextFloat())
					{
						samplestate[n] = true;
						samplecount[n]++;
					}
				}
			}
			for(int i=0; i<nodes.length;i++)
			{
				nodes[i].setMarginal(time, ((float)samplecount[i]/(float)number));
			}
		}// end for each time

	}

}
