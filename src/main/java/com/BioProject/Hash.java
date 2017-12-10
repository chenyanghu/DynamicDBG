package com.BioProject;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Set;
import java.util.TreeSet;

import org.biojava.nbio.sequencing.io.fastq.Fastq;
import org.biojava.nbio.sequencing.io.fastq.FastqReader;
import org.biojava.nbio.sequencing.io.fastq.SangerFastqReader;

import com.BioProject.MinimalPerfectHash.LongHash;

import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.HashMap;


public class Hash {
	int base;
	Set<String> kmers;
	int k;

	int P;
	MinimalPerfectHash<Long> minimalHash;
	public Hash(int k,Set<String> kmers) throws IOException{
		this.k=k;
		this.kmers=kmers;
		boolean isInjective=false;

		this.P=get_prime(k-1,kmers.size());
		//int length=kmers.size();
		int fails=0;
		Set<Long> set=new TreeSet<Long>();
		while(!isInjective){
			fails++;
			if(fails==5) {
				this.P=firstLargerPrime(this.P*2);
				fails=0;
			}
			isInjective=true;
			set.clear();
			this.base=(int)(Math.random()*(P-1));
			for(String key:kmers){
				//String key=kmers.get(i);
				long value=karpRabin(key,P);
				if(set.contains(value)){
					isInjective=false;
					//map1.clear();
					break;
				}
				else{
					set.add(value);
					//map1.put(key, value);
				}
			}
			
		}
		LongHash l = new LongHash();;
		byte[] desc=MinimalPerfectHash.generate(set, l);
		minimalHash=new MinimalPerfectHash<Long>(desc,l);
		/*String filepath="/cise/homes/jinhao/P2P/Log/log_peer_" + 1 + ".log";
		File place=new File("/cise/homes/jinhao/P2P/Log/peer_"+1+"/");
		//System.out.println(p.id+String.valueOf(place.exists()));
		if(!place.exists()){
			//System.out.println("create"+"/P2P/peer_"+p.id+"/");
			place.mkdirs();
		}
		BufferedWriter writefile=new BufferedWriter(new FileWriter(filepath,true));
		for(String key:kmers){
			 //System.out.println(key+":"+hmap2.get(hmap.get(key)));
			 writefile.write(key+":"+minimalHash.get(karpRabin(key,P))+"\n");
		}
		writefile.close();*/
		/*long count=0;
		for(String key:map1.keySet()){
			map2.put(map1.get(key), count);
			count++;
		}*/
	}
	public long getIndex(String kmer){
		/*if(!map1.containsKey(kmer)){
			return -1;
		}
		long m1=map1.get(kmer);
		if(!map2.containsKey(m1)){
			return -2;
		}
		long m2=map2.get(m1);*/
		long m2=minimalHash.get(karpRabin(kmer,P));
		return m2;
	}

    public static void main(String[] args) throws FileNotFoundException,
            IOException {

        FileInputStream inputFastq = new FileInputStream("ERR243027.filt.fastq");
        FastqReader qReader = new SangerFastqReader();

        int k = 52;
        Set<String> nodes = new HashSet<String>();
        for (Fastq fastq : qReader.read(inputFastq)) {
            String read = fastq.getSequence();
//            System.out.println(read);
            for (int i = 0; i <= read.length()-k; i++) {
                nodes.add(read.substring(i, i + k - 1));
                nodes.add(read.substring(i + 1, i + k));
            }
        }
        boolean isInjective=false;
        HashMap<String, Integer> hmap = new HashMap<String, Integer>();
		int P=get_prime(k-1,nodes.size());
		//int length=kmers.size();
		Set<Long> set=new TreeSet<Long>();
		int base1 = 0;
		while(!isInjective){
			isInjective=true;
			set.clear();
			 base1=(int)(Math.random()*(P-1));
			//System.out.println(P);
			for(String key:nodes){
				//String key=kmers.get(i);
				long value=karp(key,P,base1);
				
				if(set.contains(value)){
					isInjective=false;
					//hmap.clear();
					
					break;
				}
				else{
					set.add(value);
					//hmap.put(key, (int)value);
				}
			}
			
		}
		LongHash l = new LongHash();;
		byte[] desc=MinimalPerfectHash.generate(set, l);
		MinimalPerfectHash minimalHash=new MinimalPerfectHash<Long>(desc,l);
		String filepath="/cise/homes/jinhao/P2P/Log/log_peer_" + 1 + ".log";
		File place=new File("/cise/homes/jinhao/P2P/Log/peer_"+1+"/");
		//System.out.println(p.id+String.valueOf(place.exists()));
		if(!place.exists()){
			//System.out.println("create"+"/P2P/peer_"+p.id+"/");
			place.mkdirs();
		}
		BufferedWriter writefile=new BufferedWriter(new FileWriter(filepath));
		for(String key:nodes){
			 //System.out.println(key+":"+hmap2.get(hmap.get(key)));
			 writefile.write(key+":"+minimalHash.get(karp(key,P,base1))+"\n");
		}
		writefile.close();
        /*String alphabet = "ACGT"; 
        int d = alphabet.length();
        int q = get_prime(30, nodes.size());
        System.out.println(q);

        HashMap<String, Integer> hmap = new HashMap<String, Integer>();
        int amx=0;
        for (String node: nodes) {
        	int kk=karp_rabin(node, q, d, alphabet);
        	amx=Math.max(amx, kk);
            hmap.put(node, karp_rabin(node, q, d, alphabet));
        }
        Set<Integer> set=new TreeSet<Integer>();
        for(String node:hmap.keySet()){
        	if(set.contains(hmap.get(node))){
        		System.out.println(hmap.get(node));
        	}
        	else{
        		set.add(hmap.get(node));
        	}
        }
        HashMap<String, Integer> hmap2 = new HashMap<String, Integer>();
        int value = 0;
        for (String node: nodes) {
            hmap2.put(node, value);
            value += 1;
        }*/

        //System.out.println(hmap.size());
        //System.out.println(amx);
    }

    public static int karp_rabin(String test, int q, int d, String alphabet){
        int M = test.length();
        int res = 0;
        for (int i = 0; i < M; i++) {
            res = (d*res + alphabet.indexOf(test.charAt(i))) % q;
        }
        return res;
    }
    private long karpRabin(String key,int P){
    	int m=key.length();
    	long h=0;
    	for(int j=0;j<m;j++){
    		 h = (base * h + key.charAt(j)) % P;
    	}
    	return h;
    }
    private static long karp(String key,int P, int base){
    	int m=key.length();
    	long h=0;
    	for(int j=0;j<m;j++){
    		 h = (base * h + key.charAt(j)) % P;
    	}
    	return h;
    }
    public static boolean isPrime(int n) {
        if (n%2==0) {
           return false;
        }
        int i = 3;
        while(i*i <= n) {
            if (n%i == 0) {
                return false;
            }
            i += 2;
        }
        return true;
    }

    public static int firstLargerPrime(int R) {
        int i = R + 1;
        while(!isPrime(i)) {
            i += 1;
        }
        return i;
    }

    public static int get_prime(int k, int length) {
        int R = Math.max(4, k*length*length);
        return firstLargerPrime(R);
    }
    
}
