package config;


import edu.uci.ics.jung.graph.Graph;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;


import org.apache.commons.collections15.Factory;
public class PavlovModel<V,E> {
	
	protected Factory<Graph<V, E>> graphFactory;
	protected Factory<V> vertexFactory;
	protected Factory<E> edgeFactory;
	private Graph<V, E> graph;
	Random mRand = new Random();
	public PavlovModel(Factory<Graph<V, E>> graphFactory, Factory<V> vertexFactory,
			Factory<E> edgeFactory, int[] mass_addEdd) {
		this.graphFactory = graphFactory;
		this.vertexFactory = vertexFactory;
		this.edgeFactory = edgeFactory;
		mRand = new Random();
		
	}
	

	public Graph<V, E>  getGraph(int num, int[] mas) {
		List<V> curVertexes = new ArrayList();
		Graph graph = graphFactory.create();
		int[] m=mas;
		for (int i = 0; i < m.length; i++) {
			V n = vertexFactory.create();
			graph.addVertex(n);
			for (int j = 0; j < m[i]; j++) {
				curVertexes.add(n);
			}
		}
		V[] mass=(V[]) graph.getVertices().toArray();
		Random rand = new Random();
		do{
			V a =curVertexes.get(rand.nextInt(curVertexes.size()));
			V b =curVertexes.get(rand.nextInt(curVertexes.size()));
			if(a!=b&&!graph.isNeighbor(a, b)){
				graph.addEdge(edgeFactory.create(),a, b);
				curVertexes.remove(a);
				curVertexes.remove(b);
			}
		}while(curVertexes.size()>1);
		return graph;

	}

	}

