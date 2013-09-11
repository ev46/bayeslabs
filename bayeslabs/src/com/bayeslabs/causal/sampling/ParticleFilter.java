package com.bayeslabs.causal.sampling;

import java.util.Random;
import com.bayeslabs.causal.model.Node;



public class ParticleFilter implements Sampler {

	private Node[] nodes;
	int number;

	
	boolean[] samplestate;
	boolean[][] wsamples; // weighted samples
	int[] samplecount; // prior counts

	
	
	public ParticleFilter(Node[] nodes, int numberofsamples)
	{
		this.nodes = nodes;
		this.number = numberofsamples;
	}
	
	
	
	public void run()
	{
		if(nodes.length < 1)return;
		Random R = new Random();
		
		samplecount = new int[nodes.length];
		
		
		
		for(int sample = 0; sample<number; sample++)
		{
			samplestate = new boolean[nodes.length];
			for(int n=0; n<nodes.length;n++)
			{
				float prob = 0f;
				if(nodes[n].isRoot())
				{// has no causes, use leak, nodes are already toposorted
					
					prob = nodes[n].getConditional(0);
				}else
				{ // this is caused by something else
					// first we create a cpt index out of current sample state
					int bit = nodes[n].getCauses().length - 1; // we are starting with high order bits
					int index = 0;
					for(int cause : nodes[n].getCauses())
					{
						if(samplestate[cause])
							index|= 1<<bit;	//set the bit
						bit--;
					}
					prob = nodes[n].getConditional(index);
				}
				if(prob >= R.nextFloat())
				{
					samplestate[n] = true;
					samplecount[n]++;
				}
			}
		}//end of prior samples
		for(int i=0; i<nodes.length;i++)
		{
			nodes[i].setMarginal(0, ((float)samplecount[i]/(float)number));
		}
		// now we sample time with weighted samples with replacement
		float prob;
		float[] weights;
		for(int time=1; time<nodes[0].getSteps(); time++)
		{
			weights = new float[number]; // this will hold the weights
			wsamples = new boolean[number][nodes.length];
			
			for(int sampleindex=0; sampleindex<number; sampleindex++)
			{
				samplestate = new boolean[nodes.length];
				weights[sampleindex] = 1f;
				// each sample goes through nodes
				for(int nodeindex=0; nodeindex<nodes.length; nodeindex++)
				{
					prob = 0f;
					if(nodes[nodeindex].isRoot())
					{
						if(nodes[nodeindex].hasObservation(time))
						{
							// using non-absolute evidence and sampling will account for Jeffery's rule,
							// also, because this is a root node we do not need to compute weights (there is no P(e|x))
							prob = nodes[nodeindex].getObservationAt(time);
							
						}else{
							// sample from (t-1) = transition probability * population at Xt (NOR with Leak)
							prob = 1.0f - ((1.0f -  (nodes[nodeindex].getMarginal(time - 1) * nodes[nodeindex].getContinuation())) * (1.0f - nodes[nodeindex].getConditional(0)));  
							
						}
					}else{
						if(nodes[nodeindex].hasObservation(time))
						{
							// here we need to compute weights for its causes
							float evidence = nodes[nodeindex].getObservationAt(time);
							prob = evidence;
							// for now, keep it simple - if evidence > .49 then TRUE, most times these values will be
							// 1 or 0 +/- 0.01..
							boolean observation = false;
							if(evidence > 0.49f)observation=true;
							int bit = nodes[nodeindex].getCauses().length - 1;
							int index = 0;
							if(!observation)// if observation is false, evidence weight should be be 1-evidence
							{
								evidence = (1 - evidence);
							}
							for(int c : nodes[nodeindex].getCauses())
							{
								index |= 1<<bit;
								bit--;
								if(samplestate[c] == observation)
								{
									weights[sampleindex] = weights[sampleindex] * (nodes[nodeindex].getConditional(index) * evidence);
								}else{
									weights[sampleindex] = weights[sampleindex] * ( 1 - (nodes[nodeindex].getConditional(index) * evidence));
								}
								index = 0; //reset (we only set individual bits)
							}
							
						}else{
							// caused by 
							int bit = nodes[nodeindex].getCauses().length - 1; // we are starting with high order bits
							int index = 0;
							for(int cause : nodes[nodeindex].getCauses())
							{
								if(samplestate[cause])
									index|= 1<<bit;	//set the bit
								bit--;
							}
							
							prob = nodes[nodeindex].getConditional(index);
							
						}
					}
					// sample this prob (even for nodes with evidence, in case they are not absolute)
					if(prob >= R.nextFloat())
					{
						samplestate[nodeindex]=true;
						wsamples[sampleindex][nodeindex] = true;
					}
				}
				
			} // end for each  sample
			// reset sample count
			samplecount = new int[nodes.length];
			// weighted re-sample with replacement
			int selected = 0;
			for(int sample = 0; selected < weights.length; sample++)
			{
				if(weights[sample] > R.nextFloat())
				{
					for(int nodeindex=0; nodeindex<nodes.length;nodeindex++)
					{
						if(wsamples[sample][nodeindex])
						{
							samplecount[nodeindex]++; // collect new counts
						}
					}
					selected++;
				}
				if(sample == weights.length)sample=0; // reset index instead of doing modulus on index
			}
			// now update marginals for current time
			for(int n = 0; n < nodes.length; n++)
			{
				nodes[n].setMarginal(time, ((float)samplecount[n] / (float)number));
			}
		}// end for each time
		 
	}

}
