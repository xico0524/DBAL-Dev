package Model;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import java.util.Random;
import java.util.StringTokenizer;

public class DBAL {
	public String inputTrainFile = null; 
	public String inputTestFile = null;
	public String separator = null;
	
	public ArrayList<Quadrup> trainData = new ArrayList<Quadrup>();
	public ArrayList<Quadrup> testData = new ArrayList<Quadrup>();
	
	public int uNum = 0;
	public int sNum = 0;
	public int tNum = 0;

	public double max_Num = 0; 
	public double min_Num = 0; 
	public double nor_pNum = 0; 
	public double max_nor_Num = 0; 
	public double min_nor_Num = 0; 
	
	
	public int trainDataNum = 0; 
	public int testDataNum = 0; 
	
	public double lambda_b; 
	public double lambda_u; 
	public double lambda_s;	
	public double lambda_t;	
	
	public double[][] U; 
	public double[][] S;
	public double[][] T;
	
	public double[]	A; 
	public double[]	A_gra; 
	public double[][] A_bin; 

	
	public double[] B; 
	public double[] B_gra; 
	public double[][] B_bin; 
	
	public double[] C; 
	public double[] C_gra; 
	public double[][] C_bin; 
	
	public double [] alpha_a; 
	public double [] alpha_b; 
	public double [] alpha_c; 
	

	public double[] everyRoundRMSE;
	public double[] everyRoundMAE;
	
	public int trainingRound = 1000; 
	public int convergenceRound = 1000; 
	
	public boolean flagRMSE = false; 
	public boolean flagMAE = false; 
	
	public double minRMSE = 100; 
	public double minMAE = 100; 
	
	public int minRMSERound = 0; 
	public int minMAERound = 0;
	
	public int delayCount = 10; 

	public void initData(String inputFile,ArrayList<Quadrup> data, int T)throws IOException
	{
		
		File input = new File(inputFile); 
		BufferedReader in = new BufferedReader(new FileReader(input)); 
		String inTemp;

		while((inTemp = in.readLine()) != null ) {
			StringTokenizer st = new StringTokenizer(inTemp,separator); 
			
			String iTemp = null;
			if(st.hasMoreTokens()) 
				iTemp = st.nextToken(); 
			
			String jTemp = null;
			if(st.hasMoreTokens()) 
				jTemp = st.nextToken();
			
			String kTemp = null;
			if(st.hasMoreTokens())
				kTemp = st.nextToken();
			
			String tValueTemp = null;
			if(st.hasMoreTokens())
				tValueTemp = st.nextToken();
		
			int uID = Integer.valueOf(iTemp); 
			int sID = Integer.valueOf(jTemp);
			int tID = Integer.valueOf(kTemp);
			
			double Value = Double.valueOf(tValueTemp);

			this.uNum = (this.uNum > uID) ? this.uNum : uID;
			this.sNum = (this.sNum > sID) ? this.sNum : sID;
			this.tNum = (this.tNum > tID) ? this.tNum : tID;
			
			this.max_Num = (this.max_Num > Value) ? this.max_Num : Value;
			this.min_Num = (this.min_Num < Value) ? this.min_Num : Value;	
			
			if(T==0) {
				this.testDataNum++;
			}
			else {
				this.trainDataNum++;
			}
			
			this.nor_pNum = Math.log10(Value+1);
			
			this.max_nor_Num = (this.max_nor_Num > this.nor_pNum) ? this.max_nor_Num : this.nor_pNum;
			this.min_nor_Num = (this.min_nor_Num < this.nor_pNum) ? this.min_nor_Num : this.nor_pNum;
			
			Quadrup qtemp = new Quadrup();
			qtemp.uID = uID;
			qtemp.sID = sID;
			qtemp.tID = tID;
			
			qtemp.value = this.nor_pNum;
		
			data.add(qtemp);
		}

		in.close();
	}
	
	public DBAL(String inputTrainFile, String inputTestFile, String separator) 
	{
		this.inputTrainFile = inputTrainFile;
		this.inputTestFile = inputTestFile;
		this.separator = separator;
	}
	
	public Map<Integer, ArrayList<RTuple>> USlice = null;
	public Map<Integer, ArrayList<RTuple>> SSlice = null;
	public Map<Integer, ArrayList<RTuple>> TSlice = null;

	public void partSlice() 
	{
		USlice = new HashMap<Integer,ArrayList<RTuple>>();
		SSlice = new HashMap<Integer,ArrayList<RTuple>>();
		TSlice = new HashMap<Integer,ArrayList<RTuple>>();
		
		for (Quadrup slice1: trainData)
		{
			if(USlice.containsKey(Integer.valueOf(slice1.uID))) 
			{
				RTuple rtemp = new RTuple();
				rtemp.rowID = slice1.sID;
				rtemp.colID = slice1.tID;
				rtemp.mvalue = slice1.value;
				USlice.get(Integer.valueOf(slice1.uID)).add(rtemp);
			}else {
				ArrayList<RTuple> uSlice = new ArrayList<RTuple>();
				RTuple rtemp = new RTuple();
				rtemp.rowID = slice1.sID;
				rtemp.colID = slice1.tID;
				rtemp.mvalue = slice1.value;
				uSlice.add(rtemp);
				USlice.put(Integer.valueOf(slice1.uID),uSlice);	
			}			
			
			if(SSlice.containsKey(Integer.valueOf(slice1.sID))) 
			{
				RTuple rtemp = new RTuple();
				rtemp.rowID = slice1.uID;
				rtemp.colID = slice1.tID;
				rtemp.mvalue = slice1.value;
				SSlice.get(Integer.valueOf(slice1.sID)).add(rtemp);
				
			}else {
				ArrayList<RTuple> sSlice = new ArrayList<RTuple>();
				RTuple rtemp = new RTuple();
				rtemp.rowID = slice1.uID;
				rtemp.colID = slice1.tID;
				rtemp.mvalue = slice1.value;
				sSlice.add(rtemp);
				SSlice.put(Integer.valueOf(slice1.sID),sSlice);	
			}
	
			if(TSlice.containsKey(Integer.valueOf(slice1.tID)))
			{
				RTuple rtemp = new RTuple();
				rtemp.rowID = slice1.uID;
				rtemp.colID = slice1.sID;
				rtemp.mvalue = slice1.value;
				TSlice.get(Integer.valueOf(slice1.tID)).add(rtemp);
				
			}else {
				ArrayList<RTuple> tSlice = new ArrayList<RTuple>();
				RTuple rtemp = new RTuple();
				rtemp.rowID = slice1.uID;
				rtemp.colID = slice1.sID;
				rtemp.mvalue = slice1.value;
				tSlice.add(rtemp);
				TSlice.put(Integer.valueOf(slice1.tID),tSlice);	
			}
			
		}
	}

	public double randMax = 0.5; 
	public double randMin = 4.9E-324; 
	public void initUST(int rank) 
	{
		U = new double[this.uNum+1][rank+1]; 
		S = new double[this.sNum+1][rank+1];
		T = new double[this.tNum+1][rank+1];
		
		Random randomu = new Random();
		for(int i=1; i<=uNum;i++) 
		{
			for(int j=1; j<=rank; j++) 
			{
				U[i][j] = randMin + randomu.nextDouble() * (randMax - randMin);
			}
		} 
		
		Random randoms = new Random();
		for(int i=1; i<=sNum; i++) 
		{
			for(int j=1; j<=rank; j++) 
			{
				S[i][j] = randMin + randoms.nextDouble() * (randMax - randMin);
			}
		}
		
		Random randomt = new Random();
		for(int i=1; i<=tNum; i++)
		{
			for(int j=1; j<=rank; j++)
			{
				T[i][j] = randMin + randomt.nextDouble() * (randMax - randMin);
			}
		}
	}

	public double biasMax = randMax;
	public double biasMin = 4.9E-324;
	
	public int rating_miu = 0; 
	
	public int item_bin = 100; 
	public int item_num = 0; 
	
	public double miu = 0; 
	
	public int seta1 = 10; 
	public int seta2 = seta1;
	public int seta3 = seta1;
	
	
	public void initBias(int rank) 
	{
		A = new double[this.uNum+1];
		B = new double[this.sNum+1];
		C = new double[this.tNum+1];
		
		double miu_t = 0;
		for(Quadrup temp : trainData) 
		{
			miu_t += temp.value;  
		}
		this.miu = miu_t / trainDataNum; 
		
		this.item_num = (this.tNum / this.item_bin) + 1;
		
		Random randoma = new Random(System.currentTimeMillis());
		for(int i=1; i<=uNum;i++) 
		{
			A[i] = biasMin + randoma.nextDouble() * (biasMax - biasMin);
		}

		alpha_a = new double[this.uNum+1];
		for(int i=1; i<=uNum; i++) {
			alpha_a[i] = biasMin + randoma.nextDouble() * (biasMax - biasMin);
		} 	
		
		A_gra = new double[this.uNum+1]; 
		for(int i=1; i<=uNum; i++)
		{
			double a_gra_Real = 0;
			int a_gra_Num = 0;
			
			if(USlice.containsKey(Integer.valueOf(i))) 
			{
				ArrayList<RTuple> a_gra_Slice = new ArrayList<RTuple>(USlice.get(Integer.valueOf(i)));
				
				for(RTuple atemp : a_gra_Slice) 
				{
					a_gra_Real += atemp.mvalue; 
					a_gra_Num++;
				}
				A_gra[i] = (a_gra_Real - miu*a_gra_Num) / (a_gra_Num + seta1);
			}
		}

		A_bin = new double[uNum+1][item_bin+1];
		for(int i=1; i<=uNum; i++)
		{
			for(int i_bin=1; i_bin<=item_bin; i_bin++) 
			{
				A_bin[i][i_bin] = biasMin + randoma.nextDouble() * (biasMax - biasMin);
			}
		}

		Random randomb = new Random(System.currentTimeMillis());
		for(int j=1; j<=sNum;j++) 
		{
			B[j] = biasMin + randomb.nextDouble() * (biasMax - biasMin);	
		}

		alpha_b = new double[this.sNum+1];
		for(int j=1; j<=sNum; j++) {
			alpha_b[j] = biasMin + randomb.nextDouble() * (biasMax - biasMin);
		} 	

		B_gra = new double[this.sNum+1]; 
		for(int j=1; j<=sNum; j++)
		{
			double b_gra_Real = 0;
			int b_gra_Num = 0;
			
			if(SSlice.containsKey(Integer.valueOf(j))) 
			{
				ArrayList<RTuple> b_gra_Slice = new ArrayList<RTuple>(SSlice.get(Integer.valueOf(j)));
				
				for(RTuple btemp : b_gra_Slice) 
				{
					b_gra_Real += btemp.mvalue;
					b_gra_Num++;
				}
				B_gra[j] = (b_gra_Real - miu*b_gra_Num) / (b_gra_Num + seta2);			
			}
		}

		B_bin = new double[sNum+1][item_bin+1];
		for(int j=1; j<=sNum; j++)
		{
			for(int j_bin=1; j_bin<=item_bin; j_bin++) 
			{
				B_bin[j][j_bin] = biasMin + randomb.nextDouble() * (biasMax - biasMin);
			}
		}

		Random randomc = new Random();
		for(int k=1; k<=tNum; k++) 
		{
			C[k] = biasMin + randomc.nextDouble() * (biasMax - biasMin);
		}
		alpha_c = new double[this.tNum+1];
		for(int k=1; k<=tNum; k++) {
			alpha_c[k] = biasMin + randomc.nextDouble() * (biasMax - biasMin);
		} 
	
		C_gra = new double[this.tNum+1]; 
		for(int k=1; k<=tNum; k++)
		{
			double c_gra_Real = 0;
			int c_gra_Num = 0;
			
			if(TSlice.containsKey(Integer.valueOf(k))) 
			{
				ArrayList<RTuple> c_gra_Slice = new ArrayList<RTuple>(TSlice.get(Integer.valueOf(k)));
				
				for(RTuple ctemp : c_gra_Slice) 
				{
					c_gra_Real += ctemp.mvalue; 
					c_gra_Num++;
				}
				C_gra[k] = (c_gra_Real - miu*c_gra_Num) / (c_gra_Num + seta3);				
			}
		}

		C_bin = new double[tNum+1][item_bin+1];
		for(int k=1; k<=tNum; k++)
		{
			for(int k_bin=1; k_bin<=item_bin; k_bin++) 
			{
				C_bin[k][k_bin] = biasMin + randomc.nextDouble() * (biasMax - biasMin);
			}
		}
				
	}

	public int NP = 10; 
	public int L = 3; 
	
	public double lambda1_up = 1.0E-2; 
	public double lambda1_low = 1.0E-3; 

	public double lambda2_up = 1.0E-1; 
	public double lambda2_low = 1.0E-2; 
	
	public double eta_up = 1.0E-2; 
	public double eta_low = 1.0E-3; 

	public double F = 0.9; 
	public double CR = 0.1;
	
	public double[][] DA = new double[NP][L]; 
	public double[][] DB = new double[NP][L];
	public double[][] DC = new double[NP][L];
	
	public void train(int rank) 
	{

		Random rand = new Random(); 
		
		for(int n=0; n<NP; n++) 
		{
			for(int l=0; l<L; l++) 
			{
				if(l==0) 
				{
					this.DA[n][l] = lambda1_low + rand.nextDouble()*(lambda1_up-lambda1_low); 
				}
				if(l==1)
				{
					this.DA[n][l] = lambda2_low + rand.nextDouble()*(lambda2_up-lambda2_low); 
				}
				if(l==2) 
				{
					this.DA[n][l] = eta_low + rand.nextDouble()*(eta_up-eta_low);
				}
			}
		}

		double[][][] U_DA = new double[NP][this.uNum+1][rank+1]; 
		double[][][] S_DA = new double[NP][this.sNum+1][rank+1]; 
		double[][][] T_DA = new double[NP][this.tNum+1][rank+1]; 
		double[][] A_DA = new double[NP][this.uNum+1];
		double[][] B_DA = new double[NP][this.sNum+1]; 
		double[][] C_DA = new double[NP][this.tNum+1]; 
		double[][] alpha_a_DA = new double[NP][this.uNum+1];
		double[][] alpha_b_DA = new double[NP][this.sNum+1];
		double[][] alpha_c_DA = new double[NP][this.tNum+1];
		double[][][] A_bin_DA = new double[NP][uNum+1][item_bin+1];
		double[][][] B_bin_DA = new double[NP][sNum+1][item_bin+1];
		double[][][] C_bin_DA = new double[NP][tNum+1][item_bin+1];
		
		double[][][] U_DC = new double[NP][this.uNum+1][rank+1]; 
		double[][][] S_DC = new double[NP][this.sNum+1][rank+1]; 
		double[][][] T_DC = new double[NP][this.tNum+1][rank+1]; 
		double[][] A_DC = new double[NP][this.uNum+1]; 
		double[][] B_DC = new double[NP][this.sNum+1]; 
		double[][] C_DC = new double[NP][this.tNum+1];
		double[][] alpha_a_DC = new double[NP][this.uNum+1];
		double[][] alpha_b_DC = new double[NP][this.sNum+1];
		double[][] alpha_c_DC = new double[NP][this.tNum+1];
		double[][][] A_bin_DC = new double[NP][uNum+1][item_bin+1];
		double[][][] B_bin_DC = new double[NP][sNum+1][item_bin+1];
		double[][][] C_bin_DC = new double[NP][tNum+1][item_bin+1];

		double starttime = System.currentTimeMillis();
		
		everyRoundRMSE = new double[trainingRound+1];
		
		everyRoundMAE = new double[trainingRound+1];
		
		minRMSE = 100;	
		minMAE =100;  
		
		minRMSERound = 0; 
		minMAERound = 0; 

		for(int tr=1; tr<=trainingRound; tr++)
		{
			double starttime1 = System.currentTimeMillis();
			
			Random rand_mutation = new Random();
			for(int n=0; n<NP; n++)
			{
				int r1 = 0, r2=0, r3=0; 
				
				while(r1==n || r2==n || r3==n ||r1==r2 || r1==r3 || r2==r3) 
				{
					r1 = rand_mutation.nextInt(NP);
					r2 = rand_mutation.nextInt(NP);
					r3 = rand_mutation.nextInt(NP);
				}
				
				for(int l=0; l<L; l++) 
				{
					DB[n][l] = DA[r1][l] + F * (DA[r2][l] - DA[r3][l]);
				}
	
			}
			
			Random rand_cross = new Random();
			for(int n=0; n<NP; n++) 
			{
				int rand_l = rand_cross.nextInt(L); 
				
				for(int l=0; l<L; l++) 
				{
					double rTemp = rand_cross.nextDouble();
					if((rTemp <= CR) || l==rand_l ) 
					{
						DC[n][l] = DB[n][l];
					}
					else 
					{
						DC[n][l] = DA[n][l];
					}
				}		
			}
			
			for(int n=0; n<NP; n++)
			{
				
				for(int i=1;i<=uNum; i++) 
				{
					A_DA[n][i] = A[i];
					A_DC[n][i] = A[i];
					
					alpha_a_DA[n][i] = alpha_a[i];
					alpha_a_DC[n][i] = alpha_a[i];
					
					for(int r=1; r<=rank;r++) 
					{
						U_DA[n][i][r] = U[i][r];
						U_DC[n][i][r] = U[i][r];
					}
					
					for(int bin=1; bin<=item_bin; bin++) 
					{
						A_bin_DA[n][i][bin] = A_bin[i][bin];
						A_bin_DC[n][i][bin] = A_bin[i][bin];
					}
				}
				
				for(int j=1;j<=sNum; j++) 
				{
					B_DA[n][j] = B[j];
					B_DC[n][j] = B[j];
					
					alpha_b_DA[n][j] = alpha_b[j];
					alpha_b_DC[n][j] = alpha_b[j];
					
					for(int r=1; r<=rank;r++) 
					{
						S_DA[n][j][r] = S[j][r];
						S_DC[n][j][r] = S[j][r];
					}
					for(int bin=1; bin<=item_bin; bin++) 
					{
						B_bin_DA[n][j][bin] = B_bin[j][bin];
						B_bin_DC[n][j][bin] = B_bin[j][bin];
					}
				}
				
				for(int k=1;k<=tNum; k++) 
				{
					C_DA[n][k] = C[k];
					C_DC[n][k] = C[k];
					
					alpha_c_DA[n][k] = alpha_c[k];
					alpha_c_DC[n][k] = alpha_c[k];
					
					for(int r=1; r<=rank;r++) 
					{
						T_DA[n][k][r] = T[k][r];
						T_DC[n][k][r] = T[k][r];
					}
					for(int bin=1; bin<=item_bin; bin++) 
					{
						C_bin_DA[n][k][bin] = C_bin[k][bin];
						C_bin_DC[n][k][bin] = C_bin[k][bin];
					}
				}
		

			}
			
			for(Quadrup temp : trainData) 
			{

				double m_vir_DA [] = new double[NP]; 
				double m_vir_DC [] = new double[NP]; 
				
				double err_DA [] = new double[NP]; 
				double err_DC [] = new double[NP];
				
				int t_bin = 0; 
				double m_real = 0;
				
				m_real = temp.value;
				t_bin = (temp.tID / this.item_num) + 1;
				
				for(int n=0; n<NP; n++) 
				{
					for(int r=1; r<=rank; r++ ) 
					{
						m_vir_DA[n] += U_DA[n][temp.uID][r] * S_DA[n][temp.sID][r] * T_DA[n][temp.tID][r];
						m_vir_DC[n] += U_DC[n][temp.uID][r] * S_DC[n][temp.sID][r] * T_DC[n][temp.tID][r];
					}

					m_vir_DA[n] += A_DA[n][temp.uID] + B_DA[n][temp.sID] + C_DA[n][temp.tID];
					m_vir_DC[n] += A_DC[n][temp.uID] + B_DC[n][temp.sID] + C_DC[n][temp.tID];

					m_vir_DA[n] += alpha_a_DA[n][temp.uID]*A_gra[temp.uID] + alpha_b_DA[n][temp.sID]*B_gra[temp.sID] + alpha_c_DA[n][temp.tID]*C_gra[temp.tID];
					m_vir_DC[n] += alpha_a_DC[n][temp.uID]*A_gra[temp.uID] + alpha_b_DC[n][temp.sID]*B_gra[temp.sID] + alpha_c_DC[n][temp.tID]*C_gra[temp.tID]; 

					m_vir_DA[n] += A_bin_DA[n][temp.uID][t_bin] + B_bin_DA[n][temp.sID][t_bin] + C_bin_DA[n][temp.tID][t_bin];
					m_vir_DC[n] += A_bin_DC[n][temp.uID][t_bin] + B_bin_DC[n][temp.sID][t_bin] + C_bin_DC[n][temp.tID][t_bin];

					err_DA[n] = m_real - m_vir_DA[n];
					err_DC[n] = m_real - m_vir_DC[n];
				}
				
				for(int n=0; n<NP; n++)
				{
					for(int r=1; r<=rank;r++) 
					{
						U_DA[n][temp.uID][r] = U_DA[n][temp.uID][r]*(1-DA[n][2]*DA[n][0]) + DA[n][2]*err_DA[n]*S_DA[n][temp.sID][r]*T_DA[n][temp.tID][r];
						S_DA[n][temp.sID][r] = S_DA[n][temp.sID][r]*(1-DA[n][2]*DA[n][0]) + DA[n][2]*err_DA[n]*U_DA[n][temp.uID][r]*T_DA[n][temp.tID][r];
						T_DA[n][temp.tID][r] = T_DA[n][temp.tID][r]*(1-DA[n][2]*DA[n][0]) + DA[n][2]*err_DA[n]*U_DA[n][temp.uID][r]*S_DA[n][temp.sID][r];
						
						U_DC[n][temp.uID][r] = U_DC[n][temp.uID][r]*(1-DC[n][2]*DC[n][0]) + DC[n][2]*err_DC[n]*S_DC[n][temp.sID][r]*T_DC[n][temp.tID][r];
						S_DC[n][temp.sID][r] = S_DC[n][temp.sID][r]*(1-DC[n][2]*DC[n][0]) + DC[n][2]*err_DC[n]*U_DC[n][temp.uID][r]*T_DC[n][temp.tID][r];
						T_DC[n][temp.tID][r] = T_DC[n][temp.tID][r]*(1-DC[n][2]*DC[n][0]) + DC[n][2]*err_DC[n]*U_DC[n][temp.uID][r]*S_DC[n][temp.sID][r];
					}				

					A_DA[n][temp.uID] = A_DA[n][temp.uID] * (1-DA[n][2]*DA[n][1]) + DA[n][2]*err_DA[n];
					B_DA[n][temp.sID] = B_DA[n][temp.sID] * (1-DA[n][2]*DA[n][1]) + DA[n][2]*err_DA[n];
					C_DA[n][temp.tID] = C_DA[n][temp.tID] * (1-DA[n][2]*DA[n][1]) + DA[n][2]*err_DA[n];
					
					A_DC[n][temp.uID] = A_DC[n][temp.uID] * (1-DC[n][2]*DC[n][1]) + DC[n][2]*err_DC[n];
					B_DC[n][temp.sID] = B_DC[n][temp.sID] * (1-DC[n][2]*DC[n][1]) + DC[n][2]*err_DC[n];
					C_DC[n][temp.tID] = C_DC[n][temp.tID] * (1-DC[n][2]*DC[n][1]) + DC[n][2]*err_DC[n];
					
					alpha_a_DA[n][temp.uID] = alpha_a_DA[n][temp.uID]*(1-DA[n][2]*DA[n][1]) + DA[n][2]*err_DA[n]*A_gra[temp.uID];
					alpha_b_DA[n][temp.sID] = alpha_b_DA[n][temp.sID]*(1-DA[n][2]*DA[n][1]) + DA[n][2]*err_DA[n]*B_gra[temp.sID];
					alpha_c_DA[n][temp.tID] = alpha_c_DA[n][temp.tID]*(1-DA[n][2]*DA[n][1]) + DA[n][2]*err_DA[n]*C_gra[temp.tID];
					
					alpha_a_DC[n][temp.uID] = alpha_a_DC[n][temp.uID]*(1-DC[n][2]*DC[n][1]) + DC[n][2]*err_DC[n]*A_gra[temp.uID];
					alpha_b_DC[n][temp.sID] = alpha_b_DC[n][temp.sID]*(1-DC[n][2]*DC[n][1]) + DC[n][2]*err_DC[n]*B_gra[temp.sID];
					alpha_c_DC[n][temp.tID] = alpha_c_DC[n][temp.tID]*(1-DC[n][2]*DC[n][1]) + DC[n][2]*err_DC[n]*C_gra[temp.tID];

					A_bin_DA[n][temp.uID][t_bin] = A_bin_DA[n][temp.uID][t_bin]*(1-DA[n][2]*DA[n][1]) + DA[n][2]*err_DA[n];
					B_bin_DA[n][temp.sID][t_bin] = B_bin_DA[n][temp.sID][t_bin]*(1-DA[n][2]*DA[n][1]) + DA[n][2]*err_DA[n];
					C_bin_DA[n][temp.tID][t_bin] = C_bin_DA[n][temp.tID][t_bin]*(1-DA[n][2]*DA[n][1]) + DA[n][2]*err_DA[n];
					
					A_bin_DC[n][temp.uID][t_bin] = A_bin_DC[n][temp.uID][t_bin]*(1-DC[n][2]*DC[n][1]) + DC[n][2]*err_DC[n];
					B_bin_DC[n][temp.sID][t_bin] = B_bin_DC[n][temp.sID][t_bin]*(1-DC[n][2]*DC[n][1]) + DC[n][2]*err_DC[n];
					C_bin_DC[n][temp.tID][t_bin] = C_bin_DC[n][temp.tID][t_bin]*(1-DC[n][2]*DC[n][1]) + DC[n][2]*err_DC[n];

				}							
			}	
			
			double [] loss_term_DA = new double[NP]; 
			double [] loss_term_DC = new double[NP];

			for(int n=0; n<NP; n++)
			{
				for(Quadrup train : trainData)
				{
					int bin = 0;
					double y_hat_DA = 0;
					double y_hat_DC = 0;

					for(int r=1; r<=rank; r++) 
					{
						y_hat_DA += U_DA[n][train.uID][r]*S_DA[n][train.sID][r]*T_DA[n][train.tID-1][r];
						y_hat_DC += U_DC[n][train.uID][r]*S_DC[n][train.sID][r]*T_DC[n][train.tID-1][r];
					}

					y_hat_DA += A_DA[n][train.uID] + B_DA[n][train.sID] + C_DA[n][train.tID];
					y_hat_DC += A_DC[n][train.uID] + B_DC[n][train.sID] + C_DC[n][train.tID];

					y_hat_DA += alpha_a_DA[n][train.uID]*A_gra[train.uID] +  alpha_b_DA[n][train.sID]*B_gra[train.sID] + alpha_c_DA[n][train.tID]*C_gra[train.tID];
					y_hat_DC += alpha_a_DC[n][train.uID]*A_gra[train.uID] +  alpha_b_DC[n][train.sID]*B_gra[train.sID] + alpha_c_DC[n][train.tID]*C_gra[train.tID];

					bin = (train.tID / this.item_num) + 1;
					y_hat_DA += A_bin_DA[n][train.uID][bin] + B_bin_DA[n][train.sID][bin] + C_bin_DA[n][train.tID][bin];
					y_hat_DC += A_bin_DC[n][train.uID][bin] + B_bin_DC[n][train.sID][bin] + C_bin_DC[n][train.tID][bin];
					
					loss_term_DA[n] +=  Math.pow(train.value - y_hat_DA, 2);
					loss_term_DC[n] +=  Math.pow(train.value - y_hat_DC, 2);
				} 		
			} 

			
			double [] RMSE_DA = new double[NP];
			double [] RMSE_DC = new double[NP];	
			
			for(int n=0; n<NP; n++) 
			{
				RMSE_DA[n] = Math.sqrt(loss_term_DA[n]/trainDataNum);
				RMSE_DC[n] = Math.sqrt(loss_term_DC[n]/trainDataNum);
				
				if(RMSE_DC[n] <= RMSE_DA[n])
				{
					for(int l=0; l<L; l++) 
					{
						DA[n][l] = DC[n][l];
					}
					
				}
			}

			double min_RMSE_DA = 100; 
			double min_RMSE_DC = 100; 
			int min_RMSE_DA_NP = 0;
			int min_RMSE_DC_NP = 0;
			for(int n=0; n<NP; n++) 
			{
				if(RMSE_DA[n]<min_RMSE_DA) 
				{
					min_RMSE_DA = RMSE_DA[n]; 
					min_RMSE_DA_NP = n;
				}
				
				if(RMSE_DC[n]<min_RMSE_DC) 
				{
					min_RMSE_DC = RMSE_DC[n]; 
					min_RMSE_DC_NP = n;
				}
			}
	
			if(RMSE_DA[min_RMSE_DA_NP] < RMSE_DC[min_RMSE_DC_NP])
			{
				for(int i=1; i<=uNum;i++)
				{
					A[i] = A_DA[min_RMSE_DA_NP][i];
					alpha_a[i] = alpha_a_DA[min_RMSE_DA_NP][i];
					
					for(int r=1; r<=rank;r++) 
					{
						U[i][r] = U_DA[min_RMSE_DA_NP][i][r];
					}
					for(int bin=1; bin<=item_bin; bin++)
					{
						A_bin[i][bin] = A_bin_DA[min_RMSE_DA_NP][i][bin];
					}
				}
				
				for(int j=1; j<=sNum;j++)
				{
					B[j] = B_DA[min_RMSE_DA_NP][j];
					alpha_b[j] = alpha_b_DA[min_RMSE_DA_NP][j];
					
					for(int r=1; r<=rank;r++) 
					{
						S[j][r] = S_DA[min_RMSE_DA_NP][j][r];
					}
					for(int bin=1; bin<=item_bin; bin++)
					{
						B_bin[j][bin] = B_bin_DA[min_RMSE_DA_NP][j][bin];
					}
				}
				
				for(int k=1; k<=tNum;k++)
				{
					C[k] = C_DA[min_RMSE_DA_NP][k];
					alpha_c[k] = alpha_c_DA[min_RMSE_DA_NP][k];
					
					for(int r=1; r<=rank;r++) 
					{
						T[k][r] = T_DA[min_RMSE_DA_NP][k][r];
					}
					for(int bin=1; bin<=item_bin; bin++)
					{
						C_bin[k][bin] = C_bin_DA[min_RMSE_DA_NP][k][bin];
					}
				}
				
			}
			else 
			{
				for(int i=1; i<=uNum;i++)
				{
					A[i] = A_DC[min_RMSE_DC_NP][i];
					alpha_a[i] = alpha_a_DC[min_RMSE_DC_NP][i];
					
					for(int r=1; r<=rank;r++) 
					{
						U[i][r] = U_DC[min_RMSE_DC_NP][i][r];
					}
					for(int bin=1; bin<=item_bin; bin++)
					{
						A_bin[i][bin] = A_bin_DC[min_RMSE_DC_NP][i][bin];
					}
				}
				
				for(int j=1; j<=sNum;j++)
				{
					B[j] = B_DC[min_RMSE_DC_NP][j];
					alpha_b[j] = alpha_b_DC[min_RMSE_DC_NP][j];
					
					for(int r=1; r<=rank;r++) 
					{
						S[j][r] = S_DC[min_RMSE_DC_NP][j][r];
					}
					for(int bin=1; bin<=item_bin; bin++)
					{
						B_bin[j][bin] = B_bin_DC[min_RMSE_DC_NP][j][bin];
					}
				}
				
				for(int k=1; k<=tNum;k++)
				{
					C[k] = C_DC[min_RMSE_DC_NP][k];
					alpha_c[k] = alpha_c_DC[min_RMSE_DC_NP][k];
					
					for(int r=1; r<=rank;r++) 
					{
						T[k][r] = T_DC[min_RMSE_DC_NP][k][r];
					}
					for(int bin=1; bin<=item_bin; bin++)
					{
						C_bin[k][bin] = C_bin_DC[min_RMSE_DC_NP][k][bin];
					}
				}
			}		

			double RMSEUp = 0; 
			double MAEUp = 0;
			for(Quadrup test : testData) 
			{
				int t_bin = 0;
				double ytemp = 0;
				for(int yr=1; yr<=rank; yr++) 
				{
					ytemp += U[test.uID][yr] * S[test.sID][yr] * T[test.tID][yr];
				}
				
				ytemp += A[test.uID] + B[test.sID] + C[test.tID];
							
				ytemp += alpha_a[test.uID]*A_gra[test.uID] + alpha_b[test.sID]*B_gra[test.sID] + alpha_c[test.tID]*C_gra[test.tID];

				t_bin = (test.tID / this.item_num) + 1;
				ytemp += A_bin[test.uID][t_bin] + B_bin[test.sID][t_bin] + C_bin[test.tID][t_bin];
				
				RMSEUp += Math.pow(test.value - ytemp, 2);
				MAEUp += Math.abs(test.value - ytemp);
			}
			
			everyRoundRMSE[tr] = Math.sqrt(RMSEUp/testDataNum);
			everyRoundMAE[tr] = MAEUp / testDataNum; 
			System.out.println(+tr+everyRoundRMSE[tr]+everyRoundMAE[tr]+"\n");
			
			if((Math.abs(everyRoundRMSE[tr]-minRMSE)>=0.00001) && (everyRoundRMSE[tr]<minRMSE))
			{
				minRMSE = everyRoundRMSE[tr];
				minRMSERound = tr;
			}
			else
			{
				if((tr - minRMSERound) >= delayCount) 
				{
					flagRMSE = true;
					if(flagMAE) 
					{
						convergenceRound = tr;
						break;
					}
				}
			}
			
			
			if((Math.abs(everyRoundMAE[tr]-minMAE)>=0.00001) && (everyRoundMAE[tr] < minMAE))
			{
				minMAE = everyRoundMAE[tr];
				minMAERound = tr;
			}
			else
			{
				if((tr - minMAERound) >= delayCount) 
				{
					flagMAE = true;
					if(flagRMSE)
					{
						convergenceRound = tr;
						break; 
					}
				}
			} 
			
			System.out.println("RMSE:"+minRMSE+minRMSERound);
			System.out.println("MAE:"+minMAE+minMAERound);
			
			double endtime1 = System.currentTimeMillis();
			System.out.println((endtime1-starttime1)/1000);
			
		}
		
		double endtime = System.currentTimeMillis();
		
		System.out.println((endtime-starttime)/60000);
		System.out.println(convergenceRound);
		System.out.println("RMSE:"+minRMSE+minRMSERound);
		System.out.println("RMSE:"+minMAE+minMAERound);	
	}

	public static void main(String[] args)throws NumberFormatException,IOException,InterruptedException {

		System.out.println("************************************");
		double all_start = System.currentTimeMillis();

		DBAL svd1 = new DBAL("./Samples/train.txt", "./Samples/test.txt", ":");

		try {
			svd1.initData(svd1.inputTrainFile,svd1.trainData,1);
			svd1.initData(svd1.inputTestFile,svd1.testData , 0);
			svd1.partSlice();
			
			int r = 20;
			svd1.initUST(r);
			svd1.initBias(r);
			svd1.train(r);
		
		} catch(IOException e) {
			e.printStackTrace();		
		}

		double all_end = System.currentTimeMillis();
		System.out.println((all_end-all_start)/60000);
		System.out.println("************************************");
	}
}
