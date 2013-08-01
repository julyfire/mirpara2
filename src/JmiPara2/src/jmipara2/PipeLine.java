
package jmipara2;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;

/**
 *
 * @author weibo
 */
public class PipeLine {

    private File filename;
    private String outfile;
    private double version=3.0;
    private int window=500;
    private int step=250;
    private int start=1;
    private int distance=60;
    private double cutoff=0.8;
    private String model="overall";
    private int level=1;
    private File workingDir=new File(".");
    private String packageDir;
    private boolean append=false;
    private double progress;
    private ArrayList<String> results=new ArrayList<String>();

    private ArrayList<SimSeq> seqList;
    private ArrayList<SimSeq> segList;
    private ArrayList<PriMiRNA> priList;
    private String[] last;
    private ArrayList<HashMap> featureList;
    private ArrayList<HashMap> trueList=new ArrayList<HashMap>();

    public PipeLine() {  
        
    }

    /**
     * the main program to predict the miRNAs of a input sequence
     * @throws IOException
     */
    public void run() throws IOException{
        initialize();

        loadData(filename);

        for(SimSeq seq : seqList){

            print("Begin to extract possible pri-miRNAs from sequence "+seq.getId()+"...");
            findPrimiRNA(seq);
            seq.setSeq(null); //free space
            print("\n"+priList.size()+" putative pri-miRNAs are found\n");

            print("Begin to scan possible miRNAs from pri-miRNAs...");
            int i=0;
            trueList=new ArrayList<HashMap>();
            for(PriMiRNA pri : priList){
                findMiRNA(pri);
                priList.set(i, null); //free space
                i++;
                print(Output.decimal(i*100.0/priList.size())+"%"+Output.backspace(Output.decimal(i*100.0/priList.size())+"%"));                
            }
            print("\n"+trueList.size()+" miRNA candidates are found\n");


            reportResult(seq);

        }
        print("All have done!\nSee results at "+(workingDir.getCanonicalPath())+"\n");

    }
    
    public void run2() throws IOException{
        initialize();
        ReadFile rf=new ReadFile(filename);
        while(rf.hasSeq()){
            SimSeq seq=rf.getOneSeq();  //each seq
            setOutfileName(workingDir, filename, seq.getId());
            trueList=new ArrayList<HashMap>();
            last=new String[0];
            recordResults(seq);
            this.append=true;
            int numOfMiRNA=0;
            print("Begin to scan possible miRNAs from sequence "+seq.getId()+"...\n");
            int length=seq.getLength();
            int n=(length-window)/step+1; //num of window, num of step = n-1
            if(n<1) n=1;
            if(length-((n-1)*step+window)>19) n+=1;
            int end=0,start=this.start-1;
            for(int i=0;i<n;i++){
                
                if(start>=length) break;
                end=start+window;
                if(end>length) end=length;
                String id=seq.getName()+"_"+(start+1)+"-"+end;
                String subseq=seq.getSeq().substring(start,end);
                print("scan region "+(start+1)+"-"+end+"...");
                SimSeq frag=new SimSeq(id,subseq);  //each frag
                frag.setStart(start+1);// count from 1
                frag.setEnd(end); //count from 1
                frag.setName(seq.getId());
                
                SLScaner sl=new SLScaner();
                ArrayList<PriMiRNA> pris=sl.scanSLoopFrom(frag);
                noRepeat(pris);
                print(pris.size()+" pri-miRNA are found...");
                int before=trueList.size();
                for(PriMiRNA pri : pris){
                    findMiRNA(pri);                
                }
                int add=trueList.size()-before;
                print("generate "+add+" miRNA...");
                //output the results in time to avoid memory leak
                if(trueList.size()>100){
                    numOfMiRNA+=trueList.size();
                    recordResults(seq);
                    trueList=new ArrayList<HashMap>();
                }
                progress=Double.parseDouble(Output.decimal((i+1)*100.0/n));
//                print(Output.decimal((i+1)*100.0/n)+"%"+Output.backspace(Output.decimal((i+1)*100.0/n)+"%"));
                print(Output.decimal((i+1)*100.0/n)+"%\n");
                start+=step;
            }
            print("\n");
           
            recordResults(seq);
            numOfMiRNA+=trueList.size();
            
            print(numOfMiRNA+" miRNA candidates are found\n\n");
            
            this.append=false;

        }
        print("All have done!\nSee results at "+(workingDir.getCanonicalPath())+"\n");
    }
    
    private void noRepeat(ArrayList<PriMiRNA> pris){
        //remove repeat pri
        Iterator it=pris.iterator();
        while(it.hasNext()){
            PriMiRNA pri=(PriMiRNA)(it.next());
            for(String id:last){
                if(pri.getId().equals(id)){
                    it.remove();
                    break;
                }
            }
        }
        
        //update last array;
        int n=pris.size();
        last=new String[n];
        int i=0;
        for(PriMiRNA pri:pris){
            last[i++]=pri.getId();
        }     
    }
    
    

    /**
     * test whether given miRNAs are at given sequences
     * @throws IOException
     */
    public void testGivenMir() throws IOException{
        initialize();
        BufferedReader br=new BufferedReader(new FileReader(getFilename()));
        String line;
        int total=0;
        int hit=0;
        while((line = br.readLine()) != null){
            if(line.equals("")) continue;
            //pri_id,pri_seq,mi1_start,mi1_size,mi2_start,mi2_size,...
            String[] entry=line.split("\t");

            PriMiRNA pri=new PriMiRNA(entry[0],entry[1]);
            if(SLScaner.hasTwoLoop(pri)){
                System.out.println("W1: "+entry[0]+" is not a hairpin structure!");
                continue;
            }
            MiGenesis parser=new MiGenesis(pri);
            parser.parsePrimiRNA();

            int end5=pri.getSeq5().length();
            int start3=end5+pri.getMidBase().length()+1;
            for(int i=2;i<entry.length;i+=2){
                int miStart=Integer.parseInt(entry[i]);
                int miSize=Integer.parseInt(entry[i+1]);
                if(miStart<=start3 && miStart+miSize-1>=end5){
                        System.out.println("W2: the given miRNA spans more than half the terminal loop of "+entry[0]+" , I cannot handle at present!");
                        continue;
                    }
                parser.maturateMiRNA(miStart-1,miSize);
                HashMap feat=parser.gatherFeatures();
                //outputSV(feat);
                SVMToolKit.judge(getModel(),feat,cutoff);
                Boolean isMir=Boolean.parseBoolean(SVMToolKit.judgeResult());
                if(isMir){
                    System.out.println("Y: miRNA is detected at site "+miStart+" of "+entry[0]);
                    hit++;
                }
                else System.out.println("N: miRNA is not detected at site "+miStart+" of "+entry[0]);
                total++;
            }
        }
        br.close();
        System.out.println("cutoff="+cutoff);
        System.out.println("accuracy="+((double)hit/total));
    }

    /**
     * read a SVM format matrix file and predict it
     * @throws FileNotFoundException
     * @throws IOException
     */
    public void testProgram() throws FileNotFoundException, IOException{
        initialize();
        BufferedReader br=new BufferedReader(new FileReader(getFilename()));
        BufferedWriter bw=new BufferedWriter(new FileWriter("test.out"));
        bw.write("given_label\tpreditc_label\tsvm_value\tprobability\n");
        String line;
        while((line = br.readLine()) != null){
            double[] results=SVMToolKit.predict_test(getModel(), line);
            StringBuilder out=new StringBuilder();
            for(int i=0;i<results.length;i++)
                out.append(results[i]).append("\t");
            out.replace(out.length()-1, out.length(), "\n");
            bw.write(out.toString());
        }
        bw.close();
        br.close();
    }


    private void loadData(File fn) throws IOException{
        
        seqList=new ReadFile(fn).getSeqs();
    }

    private void findPrimiRNA(SimSeq seq){
        SLScaner sl=new SLScaner(seq);
        sl.setWindow(window);
        sl.setStep(step);
        sl.setDistance(distance);
        sl.slidingWindow();
        sl.scanStemLoop();
        sl.noRepeat();
        priList=sl.getPriList();
        sl=null; // free space
    }

    private void findMiRNA(PriMiRNA priRNA) throws IOException{ 
        MiGenesis mg=new MiGenesis(priRNA);
        try{
            mg.parsePrimiRNA();
        }catch(Exception e){
            System.out.print("Cannot parse the pri-miRNA ");
            System.out.println(priRNA.getId());
            System.out.println("sequence: "+priRNA.getSeq());
            System.out.println("structure:"+priRNA.getStr());
            return;
        }
        
        for(int w=20;w<25;w++){//miRNA size
            int i; //i refer to the start of miRNA, count from 0
           
            //5' strand
            for(i=0;i<priRNA.getTerminalLoopStart()-w;i++){//miRNA start point
                mg.maturateMiRNA(i,w);
                judgeIfMiRNA(mg.gatherFeatures());
            }
            //3' strand
            for(i=priRNA.getTerminalLoopEnd()+1;i<priRNA.getLength()-w;i++){//miRNA start point
                mg.maturateMiRNA(i,w);
                judgeIfMiRNA(mg.gatherFeatures());
            }

        }
        
        

    }



    private Boolean judgeIfMiRNA(HashMap feat) throws IOException{
        String result;        
        if(FeatureRange.rangeFilter(feat,getModel())){
            SVMToolKit.judge(getModel(),feat, cutoff);
            result=SVMToolKit.judgeResult();
            if(result.equals("TRUE"))
                trueList.add(feat);
        }
        else
           result="FALSE";
//        print(feat.get("miRNA_id")+"......"+result+Output.backspace(feat.get("miRNA_id")+"......"+result));
        return Boolean.parseBoolean(result);
    }
    
    public void setOutfileName(File dir, File infile, String seqname){
        dir.mkdirs();
        String basename=infile.getName().replaceAll("\\.\\w+", "");
//        if(dir.endsWith("/"))  dir=dir.substring(dir.length()-1);
        outfile=dir+"/"+basename+"_"+seqname;
        results.add(outfile);
        Output.fname=outfile;
    }

    /**
     * output
     * @param seq
     * @throws IOException
     */
    public void reportResult(SimSeq seq) throws IOException{
        Output.detail(trueList,seq,append);
        Output.overall(trueList,seq,append);
//        Output.distribution(trueList,seq);

    }
    
    public void recordResults(SimSeq seq) throws IOException {
        Output.detail(trueList,seq,append);
        Output.overall(trueList,seq,append);
    }

    /**
     * output a support vector
     * @param feat
     */
    public void outputSV(HashMap feat){
        double[] feats=SVMToolKit.featValueOf(model,feat);
        StringBuilder pLine=new StringBuilder();
        for(int fs=1;fs<=feats.length;fs++)
                pLine.append("\t").append(fs).append(":").append(feats[fs-1]);
        System.out.println(pLine.toString());
    }

    /**
     * should call the method first
     * @throws IOException
     */
    private void initialize() throws IOException{
        //set the work directory
        //setWorkingDir(PathSet.getWorkingDir());
        //load RNAFold library
        //in JAR package
        packageDir=PathSet.getPackageDir();
        PathSet.setLibDir(packageDir+"/lib");
        File mirbaseData=new File(PathSet.getPackageDir()+"/data/mirbase.dat");
        //in NetBeans
//        packageDir="/home/wb/Desktop/program/jmi2_v13";
//        PathSet.setLibDir("/home/wb/Desktop/program/jmi2_v13/lib");
//        File mirbaseData=new File("/home/wb/Desktop/program/jmi2_v13/data/mirbase.dat");
        
        //load SVM model
        String modelName=packageDir+"/models/"+getModel()+"_"+getLevel()+".model";
//        System.out.println(modelName);
        SVMToolKit.loadModel(modelName);

        //load miRBase data
        Output.loadMirBaseData(mirbaseData);       
    }

    public void print(String s){
        System.out.print(s);
        System.out.flush();
    }

    /**
     * @return the fname
     */
    public String getFilename() throws IOException {
        return filename.getCanonicalPath();
    }

    /**
     * @param fname the fname to set
     */
    public void setFilename(File filename) {
        this.filename = filename;
    }

    /**
     * @return the window
     */
    public int getWindow() {
        return window;
    }

    /**
     * @param window the window to set
     */
    public void setWindow(int window) {
        this.window = window;
    }

    /**
     * @return the step
     */
    public int getStep() {
        return step;
    }

    /**
     * @param step the step to set
     */
    public void setStep(int step) {
        this.step = step;
    }

    /**
     * @return the threshold
     */
    public int getDistance() {
        return distance;
    }

    /**
     * @param threshold the threshold to set
     */
    public void setDistance(int distance) {
        this.distance = distance;
    }

    /**
     * @return the model
     */
    public String getModel() {
        return model;
    }

    /**
     * @param model the model to set
     */
    public void setModel(String model) {
        this.model = model;
    }

    /**
     * @return the level
     */
    public int getLevel() {
        return level;
    }

    /**
     * @param level the level to set
     */
    public void setLevel(int level) {
        this.level = level;
    }

    /**
     * @return the workingDir
     */
    public File getWorkingDir() {
        return workingDir;
    }

    /**
     * @param workingDir the workingDir to set
     */
    public void setWorkingDir(String workingDir) {
        this.workingDir = new File(workingDir);
    }
    
    public void setWorkingDir(File workingDir){
        this.workingDir=workingDir;
    }

    /**
     * @return the cutoff
     */
    public double getCutoff() {
        return cutoff;
    }

    /**
     * @param cutoff the cutoff to set
     */
    public void setCutoff(double cutoff) {
        this.cutoff = cutoff;
    }

    /**
     * @return the progress
     */
    public double getProgress() {
        return progress;
    }

    /**
     * @param progress the progress to set
     */
    public void setProgress(double progress) {
        this.progress = progress;
    }

    /**
     * @return the results
     */
    public ArrayList<String> getResults() {
        return results;
    }

    /**
     * @param results the results to set
     */
    public void setResults(ArrayList<String> results) {
        this.results = results;
    }

}
