package com.bayeslabs.causal.model;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.TreeMap;


public class Plan {
private TreeMap<Integer, Event> events; // adjacency effect list
	
	private String key; // primary key for persistence, it's null for new models
	private String name;
	private int counter = 0;
	private Node[] nodes;
	
	
	public Plan(String name)
	{
		this.name=name;
		events = new TreeMap<Integer,Event>();
		
	}
	/**
	 * Method adds a new event (random variable) to the plan.
	 * @param name
	 */
	public void addEvent(String name)
	{
		counter++;
		events.put(counter,new Event(counter,name));
		
	}
	public void addEvent(String name, int id)
	{
		events.put(id, new Event(id,name));
	}
	/**
	 * Method used to add an existing Event to the plan.  This method is used in retrieving stored plans.
	 * @param event Event
	 */
	public void addEvent(Event event)
	{
		events.put(event.getId(), event);
	}
	
	
	/**
	 * Method connects two events with a probabilistic causal relationship.
	 * @param from specify from event id
	 * @param to specify to event id
	 */
	public void addLink(int from, int to)
	{
		events.get(to).addCause(from);
		events.get(from).addEffect(to);
	}
	/**
	 * Method sets a probability for a single or a combination of causes.
	 * @param event variable id
	 * @param index translates into a Power Set representation of active causes, e.g. 101 (cause 1 and 3 are active)
	 * @param probability causal probability
	 */
	public void setElicitation(int event, int index, float probability)
	{
		events.get(event).addElicittion(index, probability);
	}
	/**
	 * Method specifies transition probability from time X to time X+1
	 * @param event
	 * @param probability
	 */
	public void setEventContinuation(int event, float probability)
	{
		events.get(event).setContinuation(probability);
	}
	/**
	 * Method returns variable name.
	 * @return
	 */
	public String getName() {
		return name;
	}
	/**
	 * Method sets variable name.
	 * @param name
	 */
	public void setName(String name) {
		this.name = name;
	}
	/**
	 * Method removes a causal link between two variables and updates all affected distributions.
	 * @param from
	 * @param to
	 */
	public void removeLink(int from, int to)
	{
		// remove 'from' as a cause for 'to'
		events.get(to).removeCause(from);
		// remove 'to' as an effect for 'from'
		events.get(from).removeEffect(to);
	}
	/**
	 * Method removes a variable from a plan and its associated causal influences.
	 * @param id
	 */
	public void removeEvent(int id)
	{
		// first remove this id as a cause for all events
		for(int effect : events.get(id).getEffects())
		{
			events.get(effect).removeCause(id);
		}
		// now remove this id as an effect of all relevant events
		for(int cause : events.get(id).getCauses())
		{
			events.get(cause).removeEffect(id);
		}
		//finally remove it from our plan
		events.remove(id);
		
	}
	
	public Event getEvent(int id)
	{
		return events.get(id);
	}
	/**
	 * Method creates a map of events with its TOPO sort indices indexed by <Variable ID> mapped to a TopoSort index.
	 * @return HashMap<EventID,TOPOindex>
	 */
	public HashMap<Integer,Integer> indexedTopoTree()
	{
		HashMap<Integer, Integer> tree = new HashMap<Integer,Integer>();  // tree<ID,TOPO Index>
		int index = 0; // toposort index
		// create a causal adjacency list
		HashMap<Integer,ArrayList<Integer>> adjlist = new HashMap<Integer,ArrayList<Integer>>();
		LinkedList<Event> queue = new LinkedList<Event>(events.values());
		
		for(Event e : queue)
		{
			adjlist.put(e.getId(), new ArrayList<Integer>());
			for(int cause : e.getCauses())
			{
				adjlist.get(e.getId()).add(cause);
			}
		}
		
		// now build the topo sort tree
		Event e;
		while(!queue.isEmpty())
		{
			e = queue.removeFirst();
			
			if(adjlist.get(e.getId()).isEmpty())
			{
				tree.put(e.getId(), index);
				index++;
				for(int effect : events.get(e.getId()).getEffects())
				{
					adjlist.get(effect).remove(new Integer(e.getId()));// remove parentless key as cause to its effect
				}
			}else{
				queue.addLast(e); // put at the end of the queue because it's not empty yet
			}
		}
		
		return tree;
		
	}
	/**
	 * Method creates a DBN using Event objects from the plan.
	 * @return DBNode[] 
	 */
	public Node[] buildDBN(int steps)
	{
		nodes = new Node[events.size()];
		HashMap<Integer,Integer> tpi = this.indexedTopoTree();
		int[] topocauses;
		Event e;
		for(int ei : events.keySet())
		{
			e = events.get(ei);
			topocauses = new int[e.getCauses().size()];
			int iter = 0;
			for(int c : e.getCauses())
			{
				topocauses[iter] = tpi.get(c);
				iter++;
			}
			nodes[tpi.get(ei)] = new Node(e,topocauses, steps);
		}
		return nodes;
	}
	
	public String lookupEventName(int id)
	{
		return events.get(id).getName();
	}
	public String toString()
	{
		return events.toString();
	}
	
	public String getEventNames()
	{
		String names = "";
		for(Event e : events.values())
		{
			names += "["+e.getName()+"] ";
		}
		return names;
	}
	
	public TreeMap<Integer, Event> getEvents() {
		return events;
	}
	
	/**
	 * Method generates a description of a plan in the form: X causes Y, X causes Z, Z causes B, C, K (stand alone)
	 * @return String description
	 */
	public String getDescription()
	{
		StringBuffer buf = new StringBuffer();
		// first get the toposorted pairs of <ID,ORDER>
		HashMap<Integer,Integer> topopairs = this.indexedTopoTree();
		Event e;
		for(int id : topopairs.keySet())
		{
			e = events.get(id);
			if(e.getEffects().size() < 1 && e.getCauses().size() < 1)
			{
				buf.append(e.getName()+", ");
			}else{
				for(int effect : e.getEffects())
				{
					buf.append(e.getName()+" causes "+events.get(effect).getName()+", ");
				}
			}
		}
		
		return buf.toString();
	}
	
	public String getKey()
	{
		return key;
	}
	public void setKey(String key)
	{
		this.key = key;
	}
	public void setCounter(int counter)
	{
		this.counter = counter;
	}
	public int nextCounter()
	{
		return ++counter;
	}
	public int getCurrentCounter()
	{
		return counter;
	}
}
