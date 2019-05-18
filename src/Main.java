import java.io.FileWriter;
import java.io.IOException;
import java.util.stream.IntStream;

public class Main {
    static double a = 5;
    static double b = -5;
    static double С = 5;
    static double lam = 1;
    static int wX = 10;
    static int wT = 1000;

    static double[][] W = new double[wX][wT];
    static double[][]  MatrixPutTX = new double[wX][wT];
    static double h = (double) 1/wX;
    static double tau = (double) Math.pow(h,2)*(0.1);

    static double ExactSolution(double x, double t){
        return  Math.pow((С*Math.exp((-lam/(2*a))*(x+lam*t)))-(2*b)/(3*lam),-2);
    }

    static double DifferentialScheme(int i,int k){
        return W[i][k]+tau*(5*((W[i-1][k]-2*W[i][k]+W[i+1][k])/(h*h))+5*Math.pow(W[i][k],0.5)*((W[i+1][k]-W[i-1][k])/(2*h)));
    }



    public static void main(String[] args) throws IOException {

        long startSerialTime = System.nanoTime();

        for(int i=0;i<W.length;i++){
            W[i][0] = ExactSolution(i*h,0);
        }

        for(int i=0;i<W[0].length;i++){
            W[0][i] = ExactSolution(0,i*tau);
            W[W.length-1][i] = ExactSolution(1,i*tau);
        }

        for (int k = 0; k < W[0].length-1; k++){
            for (int i = 1; i < W.length-1; i++){
                W[i][k+1] = DifferentialScheme(i,k);
            }
        }
        long endSerialTime = System.nanoTime();
        long SerialTime = endSerialTime - startSerialTime;


        FileWriter filePutTX = new FileWriter("MatrixPutTX.txt");

        filePutTX.write("ListPlot3D[{"+"{"+(0.0)+","+(0.0)+","+W[0][0]+"}");
        for (int i = 0; i < W.length; i++)
        {
            for (int j = 0; j < W[0].length; j++) {
                filePutTX.write(",{" + (i * h) + "," + (j * tau) + "," + W[i][j] + "}");
            }
        }
        filePutTX.write("}, Mesh -> All]");
        filePutTX.close();

        FileWriter FileSerial = new FileWriter("serial.txt");

        double Absolute = 0;
        double Relative = 0;
        FileSerial.write("ListPlot3D[{"+"{"+(0.0)+","+(0.0)+","+W[0][0]+"}");
        for(int i=0;i<W.length;i++){
            for(int j=1;j<W[0].length;j++){
                //         file.write(String.format("%.3f",W[i][j])+"("+i+"  "+j+"), ");
                FileSerial.write(",{"+(i*h)+","+(j*tau)+","+W[i][j]+"}");
                if(Absolute < Math.abs(W[i][j]-ExactSolution(i*h,j*tau))){
                    Absolute = Math.abs(W[i][j]-ExactSolution(i*h,j*tau));
                    Relative = Absolute/W[i][j]*100;
                }
            }
        }
        FileSerial.write("}, Mesh -> All]");
        FileSerial.close();

        long startParallelTime = System.nanoTime();
        IntStream.range(0,W.length).parallel().forEach(i-> W[i][0] = ExactSolution(i*h,0));

        IntStream.range(0,W[0].length).parallel().forEach(i->{
            W[0][i] = ExactSolution(0,i*tau);
            W[W.length-1][i] = ExactSolution(1,i*tau);
        });
        for (int k = 0; k < W[0].length-1; k++) {
            int finaly = k;
            IntStream.range(1, W.length - 1).parallel().forEach((i) -> W[i][finaly + 1] = DifferentialScheme(i, finaly));
        }

        long endTime   = System.nanoTime();
        long totalTime = endTime - startParallelTime;

        FileWriter FileParallel = new FileWriter("parralel.txt");

        double AbsoluteParallel = 0;
        double RelativeParallel = 0;

        FileParallel.write("ListPlot3D[{"+"{"+(0.0)+","+(0.0)+","+W[0][0]+"}");
        for(int i=0;i<W.length;i++){
            for(int j=1;j<W[0].length;j++){
                FileParallel.write(",{"+(i*h)+","+(j*tau)+","+W[i][j]+"}");
                if(AbsoluteParallel < Math.abs(W[i][j]-ExactSolution(i*h,j*tau))){
                    AbsoluteParallel = Math.abs(W[i][j]-ExactSolution(i*h,j*tau));
                    RelativeParallel = AbsoluteParallel/W[i][j]*100;
                }
            }
        }
        FileParallel.write("}, Mesh -> All]");
        FileParallel.close();

        System.out.println("При послідовному вирішенні: ");
        System.out.println("Програма виконувалась " + SerialTime + " наносекунд");
        System.out.println("Абсолютна похибка = "+ Absolute);
        System.out.println("Відносна = "+ Relative+"\n");
        System.out.println("При паралельному вирішенні:");
        System.out.println("При паралельному вирішенні програма виконувалась " + totalTime + " наносекунд");
        System.out.println("Абсолютна похибка = "+ AbsoluteParallel);
        System.out.println("Відносна похибка  = "+ RelativeParallel);
    }
}

