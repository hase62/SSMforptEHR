package Hasegawa.IO;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

public class Writer {
	
	public static void write_AddHeader_AddRowName(String fileName, double[][] A, String[] Header, String Row[]) throws IOException {
		String tor = "\t";
		write_AddHeader_AddRowName(fileName, A, tor, Header, Row);
	}
	
    public static void write_AddHeader_AddRowName(String fileName, double[][] A, String tor, String[] Header, String[] Row) throws IOException {
        FileWriter writer = null;
        try {
            writer = new FileWriter(fileName);
            writer.write(tor);
            for (int j = 0; j < Header.length; j++) {
            	if(j!=Header.length-1) writer.write(String.format(Header[j]) + tor);
            	else writer.write(String.format(Header[j]));
            }
            writer.write(String.format("\n"));
            
            for (int i = 0; i < A.length; i++) {
            	writer.write(Row[i] + tor);
                for (int j = 0; j < A[0].length; j++){ 
                	if(j!=A[0].length-1)writer.write(String.format("%.15f", A[i][j]) + tor);
                	else writer.write(String.format("%.15f", A[i][j]));
                }
                writer.write(String.format("\n"));
                	//writer.write(A[i][j] + "\t");
                //writer.write("\n");
            }
        } catch (IOException e) {
            throw new RuntimeException("File writing error!", e);
        } finally {
            try {
                writer.close();
            } catch (IOException e) {
                throw new RuntimeException("Stream closing error!", e);
            }
        }
        writer.close();
    }
    
    public static void write_AddRowName(String fileName, double[][] A, String Row[]) throws IOException {
		String tor = "\t";
		write_AddRowName(fileName, A, tor, Row);
	}
	
    public static void write_AddRowName(String fileName, double[][] A, String tor, String[] Row) throws IOException {
        FileWriter writer = null;
        try {
            writer = new FileWriter(fileName);
            for (int i = 0; i < A.length; i++) {
            	writer.write(Row[i] + tor);
                for (int j = 0; j < A[0].length; j++){ 
                	if(j!=A[0].length-1)writer.write(String.format("%.15f", A[i][j]) + tor);
                	else writer.write(String.format("%.15f", A[i][j]));
                }
                writer.write(String.format("\n"));
                	//writer.write(A[i][j] + "\t");
                //writer.write("\n");
            }
            writer.close();
        } catch (IOException e) {
            throw new RuntimeException("File writing error!", e);
        }
    }
	
	public static void write_AddHeader(String fileName, double[][] A, String[] Header) throws IOException {
		String tor = "\t";
		write_AddHeader(fileName, A, tor, Header);
	}
	
    public static void write_AddHeader(String fileName, double[][] A, String tor, String[] Header) throws IOException {
        FileWriter writer = null;
        try {
            writer = new FileWriter(fileName);
            for (int j = 0; j < Header.length; j++) {
            	if(j!=Header.length-1) writer.write(String.format(Header[j]) + tor);
            	else writer.write(String.format(Header[j]));
            }
            writer.write(String.format("\n"));
            
            for (int i = 0; i < A.length; i++) {
                for (int j = 0; j < A[0].length; j++){ 
                	if(j!=A[0].length-1)writer.write(String.format("%.15f", A[i][j]) + tor);
                	else writer.write(String.format("%.15f", A[i][j]));
                }
                writer.write(String.format("\n"));
                	//writer.write(A[i][j] + "\t");
                //writer.write("\n");
            }
        } catch (IOException e) {
            throw new RuntimeException("File writing error!", e);
        } finally {
            try {
                writer.close();
            } catch (IOException e) {
                throw new RuntimeException("Stream closing error!", e);
            }
        }
        writer.close();
    } 
	
	public static void write(String fileName, double[][] A) throws IOException {
		write(fileName, A, "\t");
	}
	
    public static void write(String fileName, double[][] A, String tor) throws IOException {
        FileWriter writer = null;
        try {
            writer = new FileWriter(fileName);

            for (int i = 0; i < A.length; i++) {
                for (int j = 0; j < A[0].length; j++) {
                	writer.write(String.format("%.15f", A[i][j]));
                	if(j+1< A[0].length)writer.write(String.format(tor));
                }
                writer.write(String.format("\n"));
                	//writer.write(A[i][j] + "\t");
                //writer.write("\n");
            }
        } catch (IOException e) {
            throw new RuntimeException("File writing error!", e);
        } finally {
            try {
                writer.close();
            } catch (IOException e) {
                throw new RuntimeException("Stream closing error!", e);
            }
        }
        writer.close();
    } 
    
    public static void write(String fileName, String[][] S) throws IOException {
    	write(fileName, S, "\t");
    }
    
    public static void write(String fileName, String[][] S, String tor) throws IOException {
        FileWriter writer = null;
        try {
            writer = new FileWriter(fileName);

            for (int i = 0; i < S.length; i++) {
                for (int j = 0; j < S[0].length; j++) {
                	writer.write(S[i][j]);
                	if(j+1< S[0].length)writer.write(String.format(tor));
                }
                writer.write("\n");
                	//writer.write(A[i][j] + "\t");
                //writer.write("\n");
            }
        } catch (IOException e) {
            throw new RuntimeException("File writing error!", e);
        } finally {
            try {
                writer.close();
            } catch (IOException e) {
                throw new RuntimeException("Stream closing error!", e);
            }
        }
        writer.close();
    } 

    public static void write(String fileName, Double[] v) throws IOException {
        FileWriter writer = null;
        try {
            writer = new FileWriter(fileName);

            for (int i = 0; i < v.length; i++) {
                        writer.write(String.format("%.15f", v[i]));
                writer.write(String.format("\n"));
                	//writer.write(v[i] + "\t");
                //writer.write("\n");
            }
        } catch (IOException e) {
            throw new RuntimeException("File writing error!", e);
        } finally {
            try {
                writer.close();
            } catch (IOException e) {
                throw new RuntimeException("Stream closing error!", e);
            }
        }
        writer.close();
    }

    
    public static void write(String fileName, double[] v) throws IOException {
        FileWriter writer = null;
        try {
            writer = new FileWriter(fileName);

            for (int i = 0; i < v.length; i++) {
                        writer.write(String.format("%.15f", v[i]));
                writer.write(String.format("\n"));
                	//writer.write(v[i] + "\t");
                //writer.write("\n");
            }
        } catch (IOException e) {
            throw new RuntimeException("File writing error!", e);
        } finally {
            try {
                writer.close();
            } catch (IOException e) {
                throw new RuntimeException("Stream closing error!", e);
            }
        }
        writer.close();
    }
    
    public static void writeArray(String fileName, ArrayList<Integer> list) throws IOException {
	    double[] v = new double[list.size()];
        for (int i = 0; i < v.length; i++) {
			v[i] = list.get(i);
		}
        write(fileName, v);
    }

    /**
     * ClassA[] list = aList.toArray(new ClassA[aList.size()]);
     */
    public static void write(String fileName, String[] s) throws IOException {
        FileWriter writer = null;
        try {
            writer = new FileWriter(fileName);

            for (int i = 0; i < s.length; i++) {
                        writer.write(s[i]);
                writer.write("\n");
            }
        } catch (IOException e) {
            throw new RuntimeException("File writing error!", e);
        } finally {
            try {
                writer.close();
            } catch (IOException e) {
                throw new RuntimeException("Stream closing error!", e);
            }
        }
        writer.close();
    }
    
    public static void write(String fileName, double v) throws IOException {
        FileWriter writer = null;
        try {
            writer = new FileWriter(fileName);
                        writer.write(String.format("%.15f", v) + "\t");
                writer.write(String.format("\n"));
                	//writer.write(v + "\t");
                //writer.write("\n");
            
        } catch (IOException e) {
            throw new RuntimeException("File writing error!", e);
        } finally {
            try {
                writer.close();
            } catch (IOException e) {
                throw new RuntimeException("Stream closing error!", e);
            }
        }
        writer.close();
    }

	public static void write(String fileName, ArrayList<ArrayList<String>> data) throws IOException {
		// TODO Auto-generated method stub
		FileWriter writer = null;
        try {
            writer = new FileWriter(fileName);
            for (int i = 0; i < data.size(); i++) {
                for (int j = 0; j < data.get(i).size(); j++) 
                	writer.write(data.get(i).get(j) + "\t");
                writer.write("\n");
            }
        } catch (IOException e) {
            throw new RuntimeException("File writing error!", e);
        } finally {
            try {
                writer.close();
            } catch (IOException e) {
                throw new RuntimeException("Stream closing error!", e);
            }
        }
        writer.close();
		
	}
}
