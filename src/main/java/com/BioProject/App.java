package com.BioProject;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

/**
 * Hello world!
 *
 */
public class App 
{
    public static void main( String[] args ) throws IOException
    {
    	String fileName=args[0];
    	int k=Integer.parseInt(args[1]);
    	//String fileName="ERR243027.filt.fastq";
    	//int k=39;


    	long startTime = System.currentTimeMillis(); 
    	deBruijnGraph graph=new deBruijnGraph(fileName,k);
    	long endTime   = System.currentTimeMillis(); 
    	long TotalTime = endTime - startTime;       
    	System.out.println("Toal construction time: "+TotalTime+"ms");
    	String u=graph.sampleEdge.substring(0, k-1);
    	String v=graph.sampleEdge.substring(1);
    	
    	System.out.println("Query edge: "+graph.sampleEdge);
    	startTime=System.currentTimeMillis();
    	boolean flag=graph.queryEdge(u, v);
    	endTime=System.currentTimeMillis();
    	TotalTime = endTime - startTime;     
    	System.out.println("Query edge time: "+TotalTime+"ms");
    	System.out.println("Query return"+":"+flag);
    	
    	System.out.println("Remove edge: "+graph.sampleEdge);
    	startTime=System.currentTimeMillis();
    	graph.removeEdge(u, v);
    	endTime=System.currentTimeMillis();
    	TotalTime = endTime - startTime;        	
    	System.out.println("Remove edge time: "+TotalTime+"ms");
    	flag=graph.queryEdge(u, v);
    	System.out.println("Query again return"+":"+flag);
    	
    	System.out.println("Add edge: "+graph.sampleEdge);
    	startTime=System.currentTimeMillis();
    	graph.addEdge(u, v);
    	endTime=System.currentTimeMillis();
    	TotalTime = endTime - startTime;        	
    	System.out.println("Add edge time: "+TotalTime+"ms");
    	flag=graph.queryEdge(u, v);
    	System.out.println("Query again return"+":"+flag);
    	//System.out.println(graph.query("GCAAAGGTATGAACCAGAGGCGAGAGCAGT"));
    	//System.out.println(graph.query("CGCCTGCCGGAAGCCTGGCAGTAACCGTTC"));
    	/*boolean flag=false;
    	while(!flag) {
    		String[] ss=generate();
    		if(graph.query(ss[0])&&graph.query(ss[1])) {
    			flag=!graph.queryEdge(ss[0], ss[1]);
        		if(flag)System.out.println(ss[0]+" "+ss[1]);
    		}
    		
        	
    	}*/
    
    	/*String u="TTTTTTTTTTTTTTTTTTTTTAAAAATTTTTTTAATTAAACAGAAGATAAT";
    	//while(!graph.query(u))u=g();
    	String v="TTTTTTTTTTTTTTTTTTTTAAAAATTTTTTTAATTAAACAGAAGATAATT";
    	System.out.println(u);
    	System.out.println(v);
    	System.out.println(graph.queryEdge(u, v));
    	graph.removeEdge(u, v);
    	System.out.println(graph.queryEdge(u, v));
    	graph.addEdge(u, v);
    	System.out.println(graph.queryEdge(u, v));*/
    	
    }
    static String[] generate() {
    	char[] s= {'A','C','G','T'};
    	StringBuilder sb=new StringBuilder();
    	Random rand=new Random();
    	for(int i=0;i<29;i++) {
    		sb.append(s[rand.nextInt(4)]);
    	}
    	String base=sb.toString();
    	String[] res=new String[2];
    	res[0]=s[rand.nextInt(4)]+base;
    	res[1]=base+s[rand.nextInt(4)];
    	return res;
    }
    static String g() {
    	char[] s= {'A','C','G','T'};
    	StringBuilder sb=new StringBuilder();
    	Random rand=new Random();
    	for(int i=0;i<30;i++) {
    		sb.append(s[rand.nextInt(4)]);
    	}
    	return sb.toString();
    }
}
