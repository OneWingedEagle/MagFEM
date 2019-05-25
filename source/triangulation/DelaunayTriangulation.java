package triangulation;

import java.awt.Point;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.TreeSet;

import fem.Edge;
import fem.Model;
import fem.Node;
import math.Vect;
import math.util;
import triangulation.DelaunayTriangulationold.GraphEdge;


public class DelaunayTriangulation
{
	 ArrayList<Triangle> _triangles;
	 ArrayList<Edge> _edges;
	 ArrayList<Node> _vertices;

	public DelaunayTriangulation(){
		
		_triangles=new  ArrayList<Triangle>();
		_edges=new  ArrayList<Edge>();
		_vertices=new  ArrayList<Node>();
	}
	
	public static void main(String[] args) throws IOException {

		String bun=util.getFile(0);

		Model model=new Model(bun);
		
		String folder=new File(model.meshFilePath).getParentFile().getPath();
		
		ArrayList<Node> vertices = new ArrayList<Node>();

		int nNodes=model.numberOfNodes;///1000;

		for(int i=1;i<=nNodes;i++){
		//	p[i-1]=model.node[i].getCoord();
			vertices.add(model.node[i]);
		}
	
		   DelaunayTriangulation dt=new    DelaunayTriangulation();
		   
		   
		   ArrayList<Triangle> traingles=dt.triangulate(vertices);
		  				


		int nRegions=1;
		Model trianglesMesh=new Model(nRegions,traingles.size(),nNodes,"triangle");

		for(int i=1;i<=nRegions;i++){


			trianglesMesh.region[i].setMaterial("mat"+i);
			trianglesMesh.region[i].setName("reg"+i);

		}


		
		for (int i = 0; i < nNodes; i++) {
			trianglesMesh.node[i+1].setCoord(vertices.get(i).getCoord());
		}
		
		for(int i=1;i<=trianglesMesh.numberOfElements;i++){
			int[] vn=new int[3];
			  Triangle t=traingles.get(i-1);
			vn[0]=t.p1.id;
			vn[1]=t.p2.id;
			vn[2]=t.p3.id;
	
			trianglesMesh.element[i].setVertNumb(vn);
		}
		
		trianglesMesh.region[1].setFirstEl(1);
		trianglesMesh.region[1].setLastEl(trianglesMesh.numberOfElements);

		
		trianglesMesh.writeMesh(folder+"\\triangles.txt");
	}
	
	
	
	final ArrayList<Triangle> triangulate(ArrayList<Node> vertices)
	{

		for(int i = 0; i < vertices.size(); ++i)
		{
			_vertices.add(vertices.get(i));
		}
		// Store the vertices locally
		//_vertices = vertices;

		// Determinate the super triangle
		double minX = vertices.get(0).getCoord(0);
		double minY = vertices.get(0).getCoord(1);
		double maxX = minX;
		double maxY = minY;

		for(int i = 0; i < vertices.size(); ++i)
		{
			//std::cout << vertices[i].id << std::endl;

			if (vertices.get(i).getCoord(0) < minX) minX = vertices.get(i).getCoord(0);
			if (vertices.get(i).getCoord(1) < minY) minY = vertices.get(i).getCoord(1);
			if (vertices.get(i).getCoord(0) > maxX) maxX = vertices.get(i).getCoord(0);
			if (vertices.get(i).getCoord(1) > maxY) maxY = vertices.get(i).getCoord(1);
		}

		final double dx = maxX - minX;
		final double dy = maxY - minY;
		final double deltaMax = util.max(dx, dy);
		final double midx = .5*(minX + maxX);
		final double midy = .5*(minY + maxY);

		 Node p1=new Node(-1,2);
		 p1.setCoord(0,midx - 20 * deltaMax);
		 p1.setCoord(1,midy - deltaMax);
		 
		 Node p2=new Node(-1,2);
		 p2.setCoord(0,midx);
		 p2.setCoord(1, midy + 20 * deltaMax);
		 
		 Node p3=new Node(-1,2);
		 p3.setCoord(0,midx + 20 * deltaMax);
		 p3.setCoord(1,midy - deltaMax);
		 
		 p1.id=vertices.size()+1;
		 p2.id=p1.id+1;
		 p3.id=p2.id+1;
	
		//std::cout << "Super triangle " << std::endl << Triangle(p1, p2, p3) << std::endl;

		// Create a list of triangles, and add the supertriangle in it
		_triangles.add(new Triangle(p1, p2, p3));


		for(int k=0;k<vertices.size();k++)
		{
			Node p=vertices.get(k);
			
			 ArrayList<Edge> polygon = new  ArrayList<Edge>();

			for(Triangle t : _triangles)
			{

				if(t.circumCircleContains(p))
				{
					t.isBad = true;

					polygon.add(t.e1);
					polygon.add(t.e2);
					polygon.add(t.e3);
				}
				else
				{
				}
				
			}
			


		ArrayList<Triangle>  triangs=new ArrayList<Triangle>();
			for(final Triangle t : _triangles)
			{
				if(!t.isBad)
					triangs.add(t);
			}

			_triangles=triangs;

			for(int k1=0; k1<polygon.size();++k1)
			{
				Edge e1=polygon.get(k1);
				
				for(int k2=k1+1; k2<polygon.size();++k2)
				{
					Edge e2=polygon.get(k2);

					if(Triangle.almost_equal(e1, e2))
					{
						e1.common = true;
						e2.common = true;
					}
				}
			}
			
			 ArrayList<Edge> polys = new  ArrayList<Edge>();
			 
				for( Edge e : polygon)
				{
					//util.pr(e.node[0].id+"  "+e.node[1].id);
					if(!e.common) 
						polys.add(e);
				}
				
				
			polygon=polys;


			for(Edge e : polygon)
				_triangles.add(new Triangle(e.node[0], e.node[1], p));

		}
util.pr(_triangles.size());

		ArrayList<Triangle> triangs=new ArrayList<Triangle>();
			for(final Triangle t : _triangles)
			{
				if( !t.containsVertex(p1) && !t.containsVertex(p2) && !t.containsVertex(p3)){

						Vect v21=t.p2.getCoord().sub(t.p1.getCoord());
						Vect v31=t.p3.getCoord().sub(t.p1.getCoord());
						Vect cross=v21.cross(v31);
						if(cross.el[2]<0){
							triangs.add(new Triangle(t.p1,t.p3,t.p2));
						}
						else

							triangs.add(t);
				}
			}

			_triangles=triangs;



		for(final Triangle t : _triangles)
		{
			_edges.add(t.e1);
			_edges.add(t.e2);
			_edges.add(t.e3);
		}

		return _triangles;
	}



}
