package Hasegawa.IO;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.StringTokenizer;

public class SettingReader {

	public ArrayList<Double> settingReader(final String str, final String tor) throws IOException{
		final ArrayList<Double> result = new ArrayList<Double>();
		final FileInputStream fInput = new FileInputStream(str);
		final BufferedReader fRead = new BufferedReader(new InputStreamReader(fInput));
		String temp;
		while((temp = fRead.readLine()) != null){
			StringTokenizer st = new StringTokenizer(temp, tor);
			String st_temp=st.nextToken();
			if(st_temp.toUpperCase().equals("TRUE")) result.add(0.0);
			else if(st_temp.toUpperCase().equals("FALSE")) result.add(1.0);
			else result.add(Double.parseDouble(st_temp));
		}
		fRead.close();
		fInput.close();
		return result;
	}
}
