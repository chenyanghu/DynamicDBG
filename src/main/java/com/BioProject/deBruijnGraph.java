package com.BioProject;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Queue;
import java.util.Set;
import java.util.TreeSet;

import org.biojava.nbio.sequencing.io.fastq.Fastq;
import org.biojava.nbio.sequencing.io.fastq.FastqReader;
import org.biojava.nbio.sequencing.io.fastq.SangerFastqReader;



public class deBruijnGraph {
	int k;
	int alpha;
	Forest forest;
	boolean[][] IN;
	boolean[][] OUT;
	int n;//size of kmers
	Hash hash;
	int size;
	public String first;
	public String second;
	String fileName;
	public String sampleEdge;
	int readCount=0;
	
	public deBruijnGraph(String fileName,int k) throws IOException{
		Set<String> kmers=new HashSet<String>();
		Set<String> edgemers=new HashSet<String>();
		this.k=k;
		this.fileName=fileName;
		long startTime = System.currentTimeMillis(); 
		ReadFastqFile(kmers,edgemers);
		long endTime   = System.currentTimeMillis(); 
    	long TotalTime = endTime - startTime;       
    	System.out.println("Read file time: "+TotalTime+"ms");
    	System.out.println("Count of reads: "+readCount);
    	System.out.println("Size of nodes: "+kmers.size());
    	startTime = System.currentTimeMillis(); 
		hash=new Hash(k,kmers);
		endTime   = System.currentTimeMillis(); 
    	TotalTime = endTime - startTime;       
    	System.out.println("Hash functions contruction time: "+TotalTime+"ms");
		n=kmers.size();
		IN=new boolean[n][4];
		OUT=new boolean[n][4];
		this.forest=new Forest(n);
		startTime = System.currentTimeMillis(); 
		for(String edge:edgemers){
			String front=edge.substring(0,k-1);
			String back=edge.substring(1,k);
			char first=edge.charAt(0);
			char last=edge.charAt(k-1);
			
			
			
			
			//System.out.println(edge);
			//System.out.println(front);
			//System.out.println(f(front));
			//System.out.println(last);
			OUT[f(front)][charToIndex(last)]=true;
			IN[f(back)][charToIndex(first)]=true;
		}
		endTime   = System.currentTimeMillis(); 
    	TotalTime = endTime - startTime;       
    	System.out.println("IN and OUT contruction time: "+TotalTime+"ms");
    	startTime = System.currentTimeMillis(); 
		construct_forest(kmers,k*2);
		endTime   = System.currentTimeMillis(); 
    	TotalTime = endTime - startTime;       
    	System.out.println("Forest contruction time: "+TotalTime+"ms");
		
	}
	int f(String kmer){
		long res=hash.getIndex(kmer);
		return (int)res;
	}
	int charToIndex(char letter){
		if(letter=='A')return 0;
		if(letter=='C')return 1;
		if(letter=='G')return 2;
		if(letter=='T')return 3;
		return 0;
	}
	char indexTochar(int index) {
		if(index==0)return'A';
		if(index==1)return'C';
		if(index==2)return'G';
		if(index==3)return'T';
		return 'N';
	}
	void ReadFastqFile(Set<String> kmers,Set<String> edgemers) throws IOException{
		FileInputStream inputFastq = new FileInputStream(this.fileName);
        FastqReader qReader = new SangerFastqReader();

 
        String read="";
        
        for (Fastq fastq : qReader.read(inputFastq)) {
            read = fastq.getSequence();
            readCount++;
//            System.out.println(read);
            for (int i = 0; i <= read.length()-k; i++) {
            	//System.out.println(kmers.size());
                kmers.add(read.substring(i, i + k - 1));
                kmers.add(read.substring(i + 1, i + k));
            }
            for (int i = 0; i < read.length()-k; i++) {
            	sampleEdge=read.substring(i, i + k );
            	edgemers.add(read.substring(i, i + k ));
            	edgemers.add(read.substring(i + 1, i + k+1));
            }
        }
	}
	void move_kmer(Set<String> kmers,Set<String> visited,String kmer){
		kmers.remove(kmer);
		visited.add(kmer);
	}
	void store(String mer){
		forest.storeNode(f(mer), mer);
	}
	
	boolean query(String m){
		int h=f(m);
		if(h<0||h>=n){
			return false;
		}
		int hop=0;
		boolean in;
		int letter;
		if(this.forest.roots.containsKey(h)){
			if(m.equals(this.forest.roots.get(h)))return true;
		}
		while(!this.forest.isStored(h)){
			hop++;
			if(hop>3*this.alpha+10)return false;
			in=this.forest.parent_in_IN(h);
			if(in)letter=charToIndex(m.charAt(k-2));
			else letter=charToIndex(m.charAt(0));
			String parent=this.forest.getNext(h, m);
		//	System.out.println(m);
			//System.out.println(parent);
			m=parent;
			h=f(m);
			if(h>=n||h<0){
				return false;
			}
			if(in){
				if(!OUT[h][letter])return false;
			}
			else{
				if(!IN[h][letter])return false;
			}
		}
		if(m.equals(this.forest.roots.get(h))){
			return true;
		}
		return false;
	}
	boolean queryEdge(String u,String v){
		int fu=f(u);
		int fv=f(v);
		int outIndex=charToIndex(v.charAt(k-2));
		int inIndex=charToIndex(u.charAt(0));
		boolean flag=OUT[fu][outIndex];
		boolean flag2=IN[fv][inIndex];
		return OUT[fu][outIndex]&&IN[fv][inIndex];
	}
	void construct_forest(Set<String> kmers,int a){
		alpha=a;
		Set<String> visited_mers = new TreeSet<String>();
		int[] h=new int[n];
		String[] p1=new String[n];
		String[] p2=new String[n];
		String[] p=new String[n];
		while(visited_mers.size()!=n&&kmers.size()!=0){
			String root=kmers.iterator().next();
			store(root);
			move_kmer(kmers,visited_mers,root);
			int r=f(root);
			p1[r]=root;
			p2[r]=root;
			h[r]=0;
			Queue<String> queue=new LinkedList<String>();
			queue.add(root);
			List<String> neighbors=new ArrayList<String>();
			List<Boolean> inorout=new ArrayList<Boolean>();
			while(!queue.isEmpty()){
				String c=queue.poll();
				getNeighbors(c,neighbors,inorout);
				char first=c.charAt(0);
				char last=c.charAt(c.length()-1);
				for(int i=0;i<neighbors.size();i++){
					String m=neighbors.get(i);
					if(!visited_mers.contains(m)){
						queue.add(m);
						move_kmer(kmers,visited_mers,m);
						int fm=f(m);
						
						int fc=f(c);
						p[fm]=c;
						if(inorout.get(i)){
							forest.setNode(fm, false, last);
						}
						else{
							forest.setNode(fm, true, first);
						}
						h[fm]=h[fc]+1;
						if(h[fm]<=alpha){
							p1[fm]=p1[fc];
							p2[fm]=p2[fc];
						}
						if(h[fm]>alpha&&h[fm]<=2*alpha){
							store(p1[fc]);
							p1[fm]=p1[fc];
							p2[fm]=p1[fc];
						}
						if(h[fm]==2*alpha+1){
							h[fm]=0;
							p1[fm]=m;
							p2[fm]=p1[fc];
						}
					}
				}
			}
 		}
		kmers=visited_mers;
	}
	
	void getNeighbors(String c,List<String> neighbors,List<Boolean> inorout){
		neighbors.clear();
		inorout.clear();
		int fc=f(c);
		for(int i=0;i<4;i++){
			if(IN[fc][i]){
				//Letter l=new Letter(i);
				char l=indexTochar(i);
				String e=forest.pushOnFront(c, l);
				neighbors.add(e);
				inorout.add(true);
			}
			if(OUT[fc][i]){
				//Letter l=new Letter(i);
				char l=indexTochar(i);
				String e=forest.pushOnBack(c, l);
				neighbors.add(e);
				inorout.add(false);
			}
		}
		if(neighbors.contains("GGCGAATGCCCGCCAGGCGATTGTGGCGTG")) {
			System.out.println(c);
		}
	}
	public Boolean addEdge(String u,String v){
		int fu=f(u);
		int fv=f(v);
		int outIndex=charToIndex(v.charAt(k-2));
		int inIndex=charToIndex(u.charAt(0));
		if(OUT[fu][outIndex]){
			System.out.println("Edge exisits");
			return false;
		}
		OUT[fu][outIndex]=true;
		IN[fv][inIndex]=true;
		
		Tree uT=getRoot(u);
		Tree vT=getRoot(v);
		Map<String,Integer> u_heights=new HashMap<String,Integer>(); 
		Map<String,Integer> v_heights=new HashMap<String,Integer>(); 
		int uH=getTreeHeightRoot(uT.root,u_heights);
		int vH=getTreeHeightRoot(vT.root,v_heights);
		if(uT.root.equals(vT.root)){
			return true;
		}
		char uL=indexTochar(inIndex);
		char vL=indexTochar(outIndex);
		mergeTrees(u,v,uL,vL,uT,vT,uH,vH,u_heights,v_heights);
		return true;
	}
	boolean removeEdge(String u,String v){
		for(int i=0;i<k-2;i++){
			if(u.charAt(i+1)!=v.charAt(i))return false;
		}
		int fu=f(u);
		int fv=f(v);
		if(fu<0||fv<0){
			return false;
		}
		int outIndex=charToIndex(v.charAt(k-2));
		int inIndex=charToIndex(u.charAt(0));
		OUT[fu][outIndex]=false;
		IN[fv][inIndex]=false;
		if(!this.forest.isStored(fu)&&this.forest.getNext(fu, u)==v){
			removeUpdate(u,v);
		}
		else if(!this.forest.isStored(fv)&&this.forest.getNext(fv, v)==u){
			removeUpdate(v,u);
		}
		
		return true;
	}
	void removeUpdate(String child,String parent){
		int fChild=f(child);
		//int fParent=f(parent);
		this.forest.storeNode(fChild, child);
		Tree root=getRoot(parent);
		Map<String,Integer> heights=new HashMap<String,Integer>();
		int height=getTreeHeightRoot(root.root,heights);
		List<String> sorted=new ArrayList<String>();
		sorted.addAll(heights.keySet());
		if(this.alpha>height){
			removalFix(sorted,heights);
		}
		Map<String,Integer> heights2=new HashMap<String,Integer>();
		Tree croot=getRoot(child);
		int childheight=getTreeHeightRoot(croot.root,heights2);
		List<String> sorted2=new ArrayList<String>();
		sorted.addAll(heights2.keySet());
		if(this.alpha>childheight){
			removalFix(sorted2,heights2);
		}
	}
	void removalFix(List<String> kmers,Map<String,Integer>heights){
		List<Integer> hashs=new ArrayList<Integer>();
		for(int i=0;i<kmers.size();i++){
			hashs.add(f(kmers.get(i)));
		}
		for(int i=0;i<kmers.size();i++){
			int hash=hashs.get(i);
			String mer=kmers.get(i);
			for(int j=0;j<4;j++){
				if(IN[hash][j]){
					//Letter l=new Letter(j);
					char l=indexTochar(j);
					String neighbor=this.forest.pushOnFront(mer, l);
					if(!hashs.contains(f(neighbor))){
						Tree root=getRoot(neighbor);
						Tree merRoot=getRoot(mer);
						//Letter nletter=new Letter(neighbor.charAt(0));
						char nletter=neighbor.charAt(0);
						//Letter tletter=new Letter(mer.charAt(k-2));
						char tletter=mer.charAt(k-2);
						Map<String,Integer> u_heights=new HashMap<String,Integer>(); 
						Map<String,Integer> v_heights=new HashMap<String,Integer>(); 
						int uH=getTreeHeightRoot(root.root,u_heights);
						int vH=getTreeHeightRoot(merRoot.root,v_heights);
						mergeTrees(neighbor,mer,nletter,tletter,root,merRoot,uH,vH,u_heights,v_heights);
						return;
					}
				}
				if(OUT[hash][j]){
					char l=indexTochar(j);
					String neighbor=this.forest.pushOnBack(mer, l);
					if(!hashs.contains(f(neighbor))){
						Tree root=getRoot(neighbor);
						Tree merRoot=getRoot(mer);
						//Letter nletter=new Letter(neighbor.charAt(0));
						char nletter=neighbor.charAt(0);
						//Letter tletter=new Letter(mer.charAt(k-2));
						char tletter=mer.charAt(k-2);
						Map<String,Integer> u_heights=new HashMap<String,Integer>(); 
						Map<String,Integer> v_heights=new HashMap<String,Integer>(); 
						int uH=getTreeHeightRoot(root.root,u_heights);
						int vH=getTreeHeightRoot(merRoot.root,v_heights);
						mergeTrees(mer,neighbor,tletter,nletter,merRoot,root,vH,uH,v_heights,u_heights);
						return;
					}
				}
			}
		}
	}
	boolean mergeTrees(String u,String v, char uL,char vL,Tree uT,Tree vT,int uH,int vH,Map<String,Integer> u_heights,Map<String,Integer> v_heights){
		int height_u=u_heights.get(u);
		int height_v=v_heights.get(v);
		if(uT.root.equals(vT.root)){
			System.out.println("Same tree");
			return false;
		}
		if(height_u>2*this.alpha&&vH<this.alpha){
			String bp=travelUp(u,this.alpha);
			this.forest.storeNode(f(bp), bp);
			reverseEdgesToRoot(u);
			this.forest.unstoreNode(f(bp));
			this.forest.setNode(f(u), false, vL);
		}
		else if(height_v>2*this.alpha&&uH<this.alpha){
			String bp=travelUp(v,this.alpha);
			this.forest.storeNode(f(bp), bp);
			reverseEdgesToRoot(v);
			this.forest.unstoreNode(f(bp));
			this.forest.setNode(f(v), false, uL);
		}
		else if(uH<this.alpha||vH<this.alpha){
			int hops=(int) (0.5*(uH+vH+1+height_u+height_v));
			String new_root;
			if(hops<=uH+height_u){
				int h=uH+height_u-hops;
				if(h>height_u){
					h=height_u;
				}
				if(h+height_v+vH+1>3*this.alpha){
					return false;
				}
				new_root=travelUp(u,h);
				reverseEdgesToRoot(new_root);
				reverseEdgesToRoot(v);
				if(!new_root.equals(uT.root)){
					this.forest.unstoreNode(f(uT.root));
				}
				this.forest.unstoreNode(f(vT.root));
				this.forest.storeNode(f(new_root), new_root);
				this.forest.setNode(f(v), true, uL);
			}
			else{
				int h=hops-uH-height_u-1;
				if(h>height_v){
					h=height_v;
				}
				if((uH+height_u+1+hops>3*this.alpha)){
					return false;
				}
				new_root=travelUp(v,h);
				reverseEdgesToRoot(new_root);
				reverseEdgesToRoot(u);
				if(!new_root.equals(vT.root)){
					this.forest.unstoreNode(f(vT.root));
				}
				this.forest.unstoreNode(f(uT.root));
				this.forest.storeNode(f(new_root), new_root);
				this.forest.setNode(f(u), false, vL);
			}
		}
		else{
			
		}
		return true;
	}
	void reverseEdgesToRoot(String node){
		int nH=f(node);
		if(this.forest.isStored(nH))return;
		String parent=getParent(node);
		int pH=f(parent);
		String grandparent;
		int gH;
		boolean in=this.forest.parent_in_IN(nH);
		while(!this.forest.isStored(pH)){
			grandparent=getParent(parent);
			boolean parent_in=this.forest.parent_in_IN(pH);
			flipEdge(pH,nH,node,in);
			node=parent;
			parent=grandparent;
			in=parent_in;
		}
		flipEdge(pH,nH,node,in);
	}
	void flipEdge(int nH,int npH,String npk,boolean in){
		this.forest.set_parent_in_IN(nH, !in);
		char l;
		if(in){
			l=npk.charAt(k-1);
			this.forest.setLetter(nH, l);
		}
		else{
			l=npk.charAt(0);
			this.forest.setLetter(nH, l);
		}
	}
	String travelUp(String node,int hops){
		int h=f(node);
		if(this.forest.isStored(h)){
			return node;
		}
		int count=0;
		while(count<hops&&!this.forest.isStored(h)){
			String parent=getParent(node);
			h=f(parent);
			node=parent;
			count++;
		}
		return node;
	}
	int getTreeHeightRoot(String root,Map<String,Integer> heights){
		//Map<String,Integer>heights=new HashMap<String,Integer>();
		heights.put(root, 0);
		List<String> children=getChildren(root);
	
		int height=0;
		
		while(children.size()!=0){
			height++;
			List<String> cc=new ArrayList<String>();
			for(int i=0;i<children.size();i++){
				heights.put(children.get(i), height);
				List<String> cn=getChildren(children.get(i));
				for(int j=0;j<cn.size();j++){
					cc.add(cn.get(j));
				}
			}
			children=cc;
		}
		return height;
	}
	List<String> getChildren(String node){
		List<String> children=new ArrayList<String>();
		List<Boolean> inorout=new ArrayList<Boolean>();
		List<String> neighbors=new ArrayList<String>();
		getNeighbors(node,neighbors,inorout);
		int n_hash;
		for(int i=0;i<neighbors.size();i++){
			n_hash=f(neighbors.get(i));
			if(n_hash==-1)System.out.println(neighbors.get(i));
			String parent=this.forest.getNext(n_hash, neighbors.get(i));
			if(node.equals(parent)&&!this.forest.isStored(n_hash)){
				children.add(neighbors.get(i));
			}
		}
		return children;
	}
	
	class Tree{
		int height;
		String root;
		public Tree(int height,String root){
			this.height=height;
			this.root=root;
		}
	}
	Tree getRoot(String node){
		int height=0;
		int fNode=f(node);
		Tree n=new Tree(height,node);
		if(this.forest.isStored(fNode)){
			//n=new Tree(height,node);
			return n;
		}
		while(!this.forest.isStored(fNode)){
			String parent=getParent(node);
			node=parent;
			fNode=f(node);
			height++;
		}
		n.root=node;
		n.height=height;
		return n;
	}
	String getParent(String kmer){
		int f=f(kmer);
		boolean in=this.forest.parent_in_IN(f);
		int letter;
		if(in){
			letter=charToIndex(kmer.charAt(k-2));
		}
		else{
			letter=charToIndex(kmer.charAt(0));
		}
		String parent_mer=this.forest.getNext(f,kmer);
		
		return parent_mer;
	}
	public static void main(String args[]) { 
		
	}
}
