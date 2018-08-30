package config;


import edu.uci.ics.jung.graph.Graph;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;


import org.apache.commons.collections15.Factory;
public class WernickeModel<V,E> {
	
	private Graph<V, E> graph;
	Random mRand = new Random();
	public WernickeModel(Graph<V, E> graph) {
		this.graph=graph;
	}

	public Graph<V, E>  getGraph(int steps, int tries) {
		Graph g=copy(graph);
		return graph;
	}

	private Graph copy(Graph<V, E> graph2) {
		Graph g=copy(graph);
		return g;
	}

}
