package com.BioProject;

import java.io.IOException;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class Forest {
	
	

	int n;
	public Map<Integer, String> roots;
	public char[] parents;
	public boolean[] inorout;
	public boolean[] rooted;
	public Forest(int size){
		
		n=size;
		roots=new HashMap<Integer,String>(n);
		parents=new char[n];
		inorout=new boolean[n];
		rooted=new boolean[n];
	}
	int nodeIndex(int i){
		return 4*i;
	}
	void setNode(int i, boolean IN, char l){
		parents[i]=l;
		inorout[i]=IN;
		rooted[i]=false;
	}
	void setLetter(int i, char l){
		parents[i]=l;
	}
	char getLetter(int i){
		return parents[i];
	}
	
	void storeNode(int i,String kmer){
		roots.put(i, kmer);
		rooted[i]=true;
	}
	boolean isStored(int i){
		return rooted[i];
	}
	boolean parent_in_IN(int i){
		return inorout[i];
	}
	void set_parent_in_IN(int i, boolean in){
		inorout[i]=in;
	}
	String getNext(int i,String kmer){
		char l=parents[i];
		//System.out.println(l.getChar()); 
		if(inorout[i]) {
			return pushOnFront(kmer,l);
		}
		
		else{
			return pushOnBack(kmer,l);
		}
	}
	void unstoreNode(int i){
		this.roots.remove(i);
		rooted[i]=false;
		
	}
	public String pushOnFront(String kmer,char l){
		
		return l+kmer.substring(0, kmer.length()-1);
	}
	public String pushOnBack(String kmer,char l){
		return kmer.substring(1)+l;
	}
	
}
