package fem;

import math.IntVect;
import math.util;

public class EdgeSet {

	public EdgeSet(){}

	public void setEdge(Model model){

		if(model.dim==2) {setEdge2D(model); return;}

		byte[][] edgeLocalNodesHexa={{6,7},{5,4},{2,3},{1,0},{6,5},{7,4},{2,1},{3,0},{6,2},{7,3},{5,1},{4,0}};
		byte[][] edgeLocalNodesTetra={{0,1},{1,2},{2,0},{0,3},{1,3},{2,3}};
		byte[][] edgeLocalNodesPryism={{0,1},{1,2},{2,0},{3,4},{4,5},{5,3},{0,3},{1,4},{2,5}};
		byte[][] edgeLocalNodesPyramid={{0,1},{3,2},{0,3},{1,2},{0,4},{1,4},{2,4},{3,4}};
		
		byte[][] edgeLocalNodes=null;
		if(model.elCode==2)
			edgeLocalNodes=edgeLocalNodesTetra;
		else if(model.elCode==3)
			edgeLocalNodes=edgeLocalNodesPryism;
		else if(model.elCode==4)
			edgeLocalNodes=edgeLocalNodesHexa;
		else if(model.elCode==5)
			edgeLocalNodes=edgeLocalNodesPyramid;


		IntVect[][] nodeNode=new IntVect[2][model.numberOfNodes+1];

		int L=6;
		int ext=1;

		for(int i=1;i<=model.numberOfNodes;i++)
		{
			nodeNode[0][i]=new IntVect(L);
			nodeNode[1][i]=new IntVect(L);
		}

		byte[] nodeNodeIndex=new byte[model.numberOfNodes+1];
		int nEdge=0;
		int n1,n2,m;
		int numbOfEdges=2*model.numberOfNodes;
		IntVect[] edgeNodes=new IntVect[2];
		edgeNodes[0]=new IntVect(numbOfEdges);
		edgeNodes[1]=new IntVect(numbOfEdges);

		boolean[][] edgeEdgeRev=new boolean[model.numberOfElements+1][model.nElEdge];

		for(int i=1;i<=model.numberOfElements;i++){
			int[] vertNumb=model.element[i].getVertNumb();
			for(int j=0;j<model.nElEdge;j++){
				n1=vertNumb[edgeLocalNodes[j][0]];
				n2=vertNumb[edgeLocalNodes[j][1]];
				if(n2<n1) { 
					int tmp=n1; n1=n2; n2=tmp;
					edgeEdgeRev[i][j]=true;
					}
				m=util.search(nodeNode[0][n1].el,n2);

				if(m<0){
					nEdge++;
					if(nEdge==edgeNodes[0].length){
						edgeNodes[0].extend(model.numberOfNodes/2);
						edgeNodes[1].extend(model.numberOfNodes/2);
					}
					edgeNodes[0].el[nEdge]=n1;
					edgeNodes[1].el[nEdge]=n2;

					if(nodeNodeIndex[n1]==nodeNode[0][n1].length-1){
						nodeNode[0][n1].extend(ext);
						nodeNode[1][n1].extend(ext);
					}
					nodeNode[0][n1].el[nodeNodeIndex[n1]]=n2;
					nodeNode[1][n1].el[nodeNodeIndex[n1]++]=nEdge;

					model.element[i].setEdgeNumb(j,nEdge);
				}
				else
					model.element[i].setEdgeNumb(j,nodeNode[1][n1].el[m]);

			}

		}

		model.numberOfEdges=nEdge;

		model.edge=new Edge[model.numberOfEdges+1];

		for(int i=1;i<=nEdge;i++){
		//	model.edge[i]=new Edge(edgeNodes[0].el[i],edgeNodes[1].el[i]);
	
			model.edge[i]=new Edge(model.node[edgeNodes[0].el[i]],model.node[edgeNodes[1].el[i]]);
			model.edge[i].setLength(model.edgeLength(i));
		//	util.pr("edge "+i+" ("+model.edge[i].node[0].id+"¨"+model.edge[i].node[1].id+")");
		}

		for(int i=1;i<=model.numberOfElements;i++){
			model.element[i].setEdgeReverse(edgeEdgeRev[i]);
			
		}
			
		model.setMinEdge();
		model.setMaxEdge();
		
		System.gc();

	
	}

	public void setEdge2D(Model model){

		int nElVert=model.nElVert;
		int nElEdge=model.nElEdge;

		int[][] elEdge=new int[model.numberOfElements+1][model.nElEdge];
		int[][] edNode=getEdgeNodesXY(model,elEdge);

		model.nXYedges=edNode.length-1;

		int ix=0;
		int[] nodeCount=new int[model.numberOfNodes+1];

		for(int i=1;i<=model.numberOfElements;i++){
			int[] vert=model.element[i].getVertNumb();
			for(int j=0;j<nElVert;j++)	{	
				if(nodeCount[vert[j]]==0)	{	
					nodeCount[model.element[i].getVertNumb(j)]=++ix;
				}
			}
		}

		model.numberOfEdges=ix;

		model.edge=new Edge[1+model.numberOfEdges+model.nXYedges];

		ix=0;
		double pp=Math.PI*2;
		for(int i=1;i<=model.numberOfNodes;i++){
			int k=nodeCount[i];
			if(k!=0){	
				//model.edge[k]=new Edge(i,i);
				model.edge[k]=new Edge(model.node[i],model.node[i]);

				model.edge[k].setLength(1.0);
			}

		}


		for(int i=1;i<=model.numberOfElements;i++){
			int[] vertNumb=model.element[i].getVertNumb();
			for(int j=0;j<nElVert;j++){
				model.element[i].setEdgeNumb(j,nodeCount[vertNumb[j]]);
			}

		}

/*		byte b2=2;
		for(int i=1;i<=model.numberOfEdges;i++){
			model.edge[i].setDirection(b2);
		}*/

		int ned=model.numberOfEdges;


		for(int i=1;i<=model.numberOfElements;i++){
			for(int j=0;j<nElEdge;j++)
				elEdge[i][j]+=ned;
			model.element[i].edgeXYNumb=elEdge[i];
		}

		for(int i=1;i<=model.nXYedges;i++)
		{
			int n1=edNode[i][0];
			int n2=edNode[i][1];

			//model.edge[i+ned]=new Edge(n1,n2);
			model.edge[i+ned]=new Edge(model.node[n1],model.node[n2]);
			double length=model.node[n1].getCoord().sub(model.node[n2].getCoord()).norm();
			model.edge[i+ned].setLength(length);
		}


	}


	public int[][] getEdgeNodesXY(Model model,int[][] elEdge ){
		return getEdgeNodes(model, 1,model.numberOfRegions,elEdge);
	}

	public int[][] getEdgeNodes(Model model,int ir1, int ir2,int[][] elEdge ){
		byte[][] arr0={{0,1},{1,2},{0,2}};
		byte[][] arr1={{1,2},{0,3},{2,3},{1,0}};
		byte[][] arr4={{6,7},{5,4},{2,3},{1,0},{6,5},{7,4},{2,1},{3,0},{6,2},{7,3},{5,1},{4,0}};
		byte[][] edgeLocalNodes=new byte[1][1];

		if(model.elCode==0){
			edgeLocalNodes=arr0;;
		}
		else if(model.elCode==1){
			edgeLocalNodes=arr1;
		}

		else if(model.elCode==4){
			edgeLocalNodes=arr4;
		}

		boolean b=true;
		if(elEdge.length==0) b=false;

		IntVect[][] nodeNode=new IntVect[2][model.numberOfNodes+1];

		int L=6;
		int ext=4;

		for(int i=1;i<=model.numberOfNodes;i++)
		{
			nodeNode[0][i]=new IntVect(L);
			nodeNode[1][i]=new IntVect(L);
		}

		byte[] nodeNodeIndex=new byte[model.numberOfNodes+1];
		int nEdge=0;
		int n1,n2,m;
		int numbOfEdges=2*model.numberOfNodes;
		IntVect[] edgeNodes=new IntVect[2];
		edgeNodes[0]=new IntVect(numbOfEdges);
		edgeNodes[1]=new IntVect(numbOfEdges);

		for(int ir=ir1;ir<=ir2;ir++)
			for(int i=model.region[ir].getFirstEl();i<=model.region[ir].getLastEl();i++){
				int[] vertNumb=model.element[i].getVertNumb();
				for(int j=0;j<model.nElEdge;j++){
					n1=vertNumb[edgeLocalNodes[j][0]];
					n2=vertNumb[edgeLocalNodes[j][1]];
					if(n2<n1) { int tmp=n1; n1=n2; n2=tmp;}
					m=util.search(nodeNode[0][n1].el,n2);

					if(m<0){
						nEdge++;
						if(nEdge==edgeNodes[0].length){
							edgeNodes[0].extend(model.numberOfNodes/2);
							edgeNodes[1].extend(model.numberOfNodes/2);
						}
						edgeNodes[0].el[nEdge]=n1;
						edgeNodes[1].el[nEdge]=n2;

						if(nodeNodeIndex[n1]==nodeNode[0][n1].length-1){
							nodeNode[0][n1].extend(ext);
							nodeNode[1][n1].extend(ext);
						}
						nodeNode[0][n1].el[nodeNodeIndex[n1]]=n2;
						nodeNode[1][n1].el[nodeNodeIndex[n1]++]=nEdge;
						if(b)
							elEdge[i][j]=nEdge;						
					}
					else{
						if(b)
							elEdge[i][j]=nodeNode[1][n1].el[m];
					}

				}

			}

		int[][] edgeNodes2=new int[nEdge+1][2];
		double min=1e40;
		double max=0;
		for(int i=1;i<=nEdge;i++){
			edgeNodes2[i][0]=edgeNodes[0].el[i];
			edgeNodes2[i][1]=edgeNodes[1].el[i];
			double length=model.node[edgeNodes2[i][0]].getCoord().sub(model.node[edgeNodes2[i][1]].getCoord()).norm();

			if(max<length) max=length;
			if(min>length) min=length;
		}

		model.minEdgeLength=min;
		model.maxEdgeLength=max;
		return edgeNodes2;
	}



}
