package math;

public class TriangCount {
	
	public static void main(){
		int[][] comb={{0,1,2},{0,2,1},{1,0,2},{1,2,0},{2,0,1},{2,1,0}};
		for(int p=1;p<2;p++){
			
			int[][] vv=new int[1000000][3];
		int sum=20;
		int nmax=sum/2;
		
		int ppx=0;
		for(int i=1;i<=nmax;i++)
			for(int j=1;j<=nmax;j++)
				for(int k=1;k<=nmax;k++){
					if(i+j+k==sum && i+j>k && i+k>j&& j+k>i){
						
						vv[ppx][0]=i;
						vv[ppx][1]=j;
						vv[ppx][2]=k;
			
						boolean distinct=true;
						for(int t=0;t<ppx;t++){
							
							for(int h=0;h<comb.length;h++){
								int i1=comb[h][0];
								int j1=comb[h][1];
								int k1=comb[h][2];
								if(vv[t][i1]==vv[ppx][0]&& vv[t][j1]==
										vv[ppx][1]&& vv[t][k1]==vv[ppx][2]){
									distinct=false;
								break;
								}
							}
							if(!distinct){
								break;
							}
						}
						
						if(distinct){
						
							ppx++;
						}
		
				}
				}

		util.pr(sum+"  --numb of tri  "+ppx);
		for(int i=0;i<-ppx;i++)
util.hshow(vv[i]);
		}
		
	}

}
