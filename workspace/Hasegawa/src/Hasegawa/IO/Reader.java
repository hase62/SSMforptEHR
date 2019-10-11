package Hasegawa.IO;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.List;
import java.util.StringTokenizer;

public class Reader {
	
	public static int CountRow(String str) throws IOException{
		final FileInputStream fInput = new FileInputStream(str);
		final BufferedReader fRead = new BufferedReader(new InputStreamReader(fInput));
		int count = 0;
		while(fRead.readLine() != null){
			count++;
		}
		fInput.close();
		return count;
	}
	
	public static double[][] toArrayDouble(ArrayList<ArrayList<Double>> l){
		double[][] data = new double[l.size()][l.get(0).size()];
		for (int i = 0; i < data.length; i++) 
			for (int j = 0; j < data[0].length; j++) 
				data[i][j] = l.get(i).get(j);
		return data;
	}
	
	public static String[][] toArrayString(ArrayList<ArrayList<String>> l){
		String[][] data = new String[l.size()][l.get(0).size()];
		for (int i = 0; i < data.length; i++) 
			for (int j = 0; j < data[0].length; j++) 
				data[i][j] = l.get(i).get(j);
		return data;
	}
	
	public static double[] toArrayDouble(List<Double> l){
		double[] data = new double[l.size()];
		for (int i = 0; i < data.length; i++) 
			data[i] = l.get(i);
		return data;
	}
	
	public static String[] toArrayString(List<String> l){
		String[] data = new String[l.size()];
		for (int i = 0; i < data.length; i++) 
			data[i] = l.get(i);
		return data;
	}
	
	public static double[][] ReadMatrixDouble(String str, String tor) throws NumberFormatException, IOException{
		final FileInputStream fInput = new FileInputStream(str);
		final BufferedReader fRead = new BufferedReader(new InputStreamReader(fInput));
		String temp;
		ArrayList<ArrayList<Double>> l = new ArrayList<ArrayList<Double>>();
		int count = 0;
		while((temp = fRead.readLine()) != null){
			StringTokenizer st = new StringTokenizer(temp, tor);
			l.add(new ArrayList<Double>());
			while(st.hasMoreTokens()) 
				l.get(count).add(Double.parseDouble(st.nextToken()));		
			count++;
		}		
		fInput.close();
		return toArrayDouble(l);
	}
	
	public static double[] ReadVectorDouble(String str, String tor) throws NumberFormatException, IOException{
		final FileInputStream fInput = new FileInputStream(str);
		final BufferedReader fRead = new BufferedReader(new InputStreamReader(fInput));
		String temp;
		List<Double> l = new ArrayList<Double>();
		while((temp = fRead.readLine()) != null){
			StringTokenizer st = new StringTokenizer(temp, tor);
			while(st.hasMoreTokens()) 
				l.add(Double.parseDouble(st.nextToken()));		
		}		
		fInput.close();
		return toArrayDouble(l);
	}
	
	public static double ReadValueDouble(String str, String tor) throws NumberFormatException, IOException{
		final FileInputStream fInput = new FileInputStream(str);
		final BufferedReader fRead = new BufferedReader(new InputStreamReader(fInput));
		double data = 0;
		String temp = fRead.readLine();
		data = Double.parseDouble(temp);
		fInput.close();
		return data;
	}
	
	public static String[] ReadVectorString(String str) throws NumberFormatException, IOException{
		ArrayList<String> data = new ArrayList<String>();
		final FileInputStream fInput = new FileInputStream(str);
		final BufferedReader fRead = new BufferedReader(new InputStreamReader(fInput));
		String temp;
		while((temp = fRead.readLine()) != null){
			data.add(temp);
		}
		fInput.close();
		return toArrayString(data);
	}
	
	public static ArrayList<String> ReadVectorStringList(String str) throws NumberFormatException, IOException{
		ArrayList<String> data = new ArrayList<String>();
		final FileInputStream fInput = new FileInputStream(str);
		final BufferedReader fRead = new BufferedReader(new InputStreamReader(fInput));
		String temp;
		while((temp = fRead.readLine()) != null){
			data.add(temp);
		}
		fInput.close();
		return data;
	}
	
	public static ArrayList<ArrayList<String>> ReadMatrixStringAsColVector(String str, String tor) throws NumberFormatException, IOException{
		final FileInputStream fInput = new FileInputStream(str);
		final BufferedReader fRead = new BufferedReader(new InputStreamReader(fInput));
		ArrayList<ArrayList<String>> data = new ArrayList<ArrayList<String>>();
		ArrayList<String> data_list = new ArrayList<String>();
		int col = 0;
		String temp;
		while((temp = fRead.readLine()) != null){
			StringTokenizer st = new StringTokenizer(temp, tor);
			col = st.countTokens();
			while(st.hasMoreTokens()){
				data_list.add(st.nextToken());
			}
		}
		fInput.close();
		for (int i = 0; i < col; i++) {
			ArrayList<String> c = new ArrayList<String>();
			for (int j = 0; j < data_list.size()/col; j++) {
				c.add(data_list.get(j * col + i));
			}
			data.add(c);
		}
		return data;
	}

	public static String[][] ReadMatrixString(String str, String tor) throws IOException {
		final FileInputStream fInput = new FileInputStream(str);
		final BufferedReader fRead = new BufferedReader(new InputStreamReader(fInput));
		String temp;
		ArrayList<ArrayList<String>> l = new ArrayList<ArrayList<String>>();
		int count = 0;
		while((temp = fRead.readLine()) != null){
			StringTokenizer st = new StringTokenizer(temp, tor);
			l.add(new ArrayList<String>());
			while(st.hasMoreTokens()) 
				l.get(count).add(st.nextToken());
			count++;
		}		
		fInput.close();
		return toArrayString(l);
	}
}