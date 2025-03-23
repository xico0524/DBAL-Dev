package Model;

import java.io.IOException;

public class standard_deviation
{

	public static void main(String[] args)throws NumberFormatException,IOException,InterruptedException {    

		System.out.println("************************************");
			
		double [] arr1 = {1,2,3,4,5,6,7,8,9,10};
		double mean_num = 0; 
		for(int i=0; i<arr1.length;i++) 
		{
			mean_num += arr1[i];
		}
		mean_num = (mean_num / arr1.length);
		System.out.println(mean_num);
		
		double sta_dev = 0; 
		for(int j=0; j<arr1.length; j++) 
		{
			sta_dev += Math.pow((arr1[j]-mean_num), 2);
		}
		sta_dev = (sta_dev / arr1.length);
		sta_dev = Math.sqrt(sta_dev);
		System.out.println(sta_dev);
		
		System.out.println("************************************");

	}
}
