/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package jmipara2;

import java.io.File;
import java.io.IOException;

/**
 *
 * @author weibo
 */
public class JmiCMD {

    private static int test=0;
    private static double version=4.0;

    public static void main(String[] args) throws IOException {
        long start = System.currentTimeMillis();
        
        printLogo();
        
        PipeLine pl=new PipeLine();

        readArguments(args,pl);

        if(test==1){
            pl.testProgram();
            System.exit(0);
        }
        if(test==2){
            pl.testGivenMir();
            System.exit(0);
        }
        
        pl.run2();

        long end = System.currentTimeMillis();
        long time = end - start;
        System.out.println("total time:" + time/1000+"s");
    }

    public static void exit_with_help(){
                printUsage();
		System.exit(0);
    }
    
    public static void printLogo(){
        System.out.println("\n    =======================================================================\n"
                + "    |    JmiPara "+version+" can detact the possible miRNAs on given sequences    |\n"
                + "    =======================================================================\n");
        System.out.println("*** Please report bugs to helloweibo@gmail.com \n");
    }
    
    public static void printUsage(){
        System.out.println("Usage: \n"
                + "[1]start GUI: java -jar JmiPara2.jar or double click the package JmiPara2\n"
                + "[2]start CMD: java -cp JmiPara2.jar jmipara2.JmiCMD [options] input_file [output_directory]\n"
                + "[3]start MirViewer: java -cp JmiPara2.jar jmipara2.MirViewer\n\n"
		+"CMD options:\n"
		+"\t-m model:     the model used in prediction, can be animal, plant, \n"
                +"\t              virus or overall (default overall)\n"
		+"\t-l level:     the value of negative_data/positive_data, 1~20(default 1)\n"
                +"\t-s step size: step size in spliting long sequences\n"
                +"\t-w window:    window size in spliting long sequences\n"
                +"\t-d distance:  the minimal length of a hairpin\n"
                +"\t-c cutoff:    the minimal possibility that testing seq is miRNA\n"
                +"\t-t test:      1: test a given svm matrix file\n"
                +"\t              2: given a hairpin sequence, test whether a given site \n"
                +"\t                 have miRNA with given size (default 0)\n"
                +"\t-h help\n\n"+
                 "input_file:        the input file, should be fasta format at this version\n"
                +"output_directory:  the output directory, default is current directory\n\n"
                +"[References]\n"
                + "MiRPara: a SVM-based software tool for prediction of most probable microRNA"
                + " coding regions in genome scale sequences. Wu Y., Wei B., Liu H., Li T., "
                + "Rayner S. BMC Bioinformatics. 2011 Apr 19; 12(1):107\n\n");
    }

    private static void exit_with_error(){
        System.out.println("use -h to see detail help");
        System.exit(1);
    }

    private static void readArguments(String[] argv, PipeLine p){
        int i;
        if(argv.length==0) exit_with_help();
        for(i=0;i<argv.length;i++){
            if (argv[i].charAt(0) != '-')	break;
            if(++i>=argv.length && !(argv[i-1].equals("-h"))){
                System.out.println("arguments error!");
                exit_with_error();
            }
            switch(argv[i-1].charAt(1)){
                case 'm': p.setModel(argv[i]);	break;
		case 'l': p.setLevel(Integer.parseInt(argv[i]));	break;
                case 's': p.setStep(Integer.parseInt(argv[i]));  break;
                case 'w': p.setWindow(Integer.parseInt(argv[i])); break;
                case 'd': p.setDistance(Integer.parseInt(argv[i])); break;
                case 'c': p.setCutoff(Double.parseDouble(argv[i])); break;
                case 't': test=Integer.parseInt(argv[i]); break;
                case 'h': exit_with_help();
		default:
			System.err.println("unknown option");
			exit_with_error();
            }
	}
        if(!(p.getModel().equals("animal") || p.getModel().equals("plant")
                || p.getModel().equals("virus") || p.getModel().equals("overall"))){
            System.err.println("unknow taxonomy!");
            exit_with_error();
        }
        if(p.getLevel()<0){
            System.err.println("level should be greater than 0");
            exit_with_error();
        }
        if(p.getWindow()<0){
            System.err.println("window should not be less than 0");
            exit_with_error();
        }
        if(p.getStep()<0){
            System.err.println("step should not be less than 0");
            exit_with_error();
        }
        if(p.getDistance()<0){
            System.err.println("distance should not be less than 0");
            exit_with_error();
        }
        if(p.getCutoff()<0 || p.getCutoff()>1){
            System.err.println("cutoff should be a value between 0 to 1");
            exit_with_error();
        }
        // determine filenames, at least two filenames
        if(i>=argv.length){
            System.err.println("lack input file!");
            exit_with_error();
        }
        p.setFilename(new File(argv[i]));
        if(i<argv.length-1)
            p.setWorkingDir(argv[i+1]);
    }

}
