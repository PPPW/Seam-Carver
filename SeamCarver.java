import java.awt.Color;

public class SeamCarver {
    // this costs 10*N^2 + 8*N memory
    // the reference only costs 4*N^2
    private int[][] temp;
    private int[][] energy;
    private double[] distTo;
    private short[][] pixTo;
    private int wc, hc;

    // create a seam carver object based on the given picture
    public SeamCarver(Picture picture) {
        // x: horizontal right; y: vertical down
        // (x, y) -> temp[y][x]
        // to make it clear, use i to iterate x, j to y.
        wc = picture.width();
        hc = picture.height();
        temp = new int[hc][wc];
        for (int i = 0; i < width(); i++) {
            for (int j = 0; j < height(); j++) {
                temp[j][i] = picture.get(i, j).getRGB();                
            }
        }
        
        // store energy into an array
        energy = new int[width()][height()];
        for (int i = 0; i < width(); i++) {
            for (int j = 0; j < height(); j++) {
                energy[i][j] = (int) energy(i, j);
            }
        }
        
    }            

    // current picture
    public Picture picture() {
        Picture pic = new Picture(width(), height());
        for (int i = 0; i < width(); i++) {
            for (int j = 0; j < height(); j++) {
                Color c = new Color(toR(temp[j][i]), toG(temp[j][i]), toB(temp[j][i]));
                pic.set(i, j, c);               
            }
        }
        return pic;
    }
                           
    // width of current picture, columns of temp
    public int width() {
        return wc;
    }   

    // height of current picture, rows of temp
    public int height() {
        return hc;
    }  

    // energy of pixel at column x and row y
    public double energy(int x, int y) {
        if (x < 0 || x >= width())  throw new IndexOutOfBoundsException("x = " + x + " is not between 0 and " + (width() - 1));
        if (y < 0 || y >= height()) throw new IndexOutOfBoundsException("y = " + y + " is not between 0 and " + (height() - 1));
        
        if (x == 0 || y == 0 || x == width() - 1 || y == height() - 1) { return 195075; }

        double delx, dely;
        int p1, p2, p3, p4;               

        // (x, y) -> temp[y][x]
        p1 = temp[y][x - 1];
        p2 = temp[y][x + 1];  
        p3 = temp[y - 1][x];
        p4 = temp[y + 1][x];      
        //StdOut.println(toR(p1));

        delx = (toR(p2) - toR(p1)) * (toR(p2) - toR(p1)) 
             + (toG(p2) - toG(p1)) * (toG(p2) - toG(p1)) 
             + (toB(p2) - toB(p1)) * (toB(p2) - toB(p1)); 
        dely = (toR(p4) - toR(p3)) * (toR(p4) - toR(p3)) 
             + (toG(p4) - toG(p3)) * (toG(p4) - toG(p3)) 
             + (toB(p4) - toB(p3)) * (toB(p4) - toB(p3)); 

        //StdOut.println(delx + dely);
        return delx + dely;
    } 

    // sequence of indices for horizontal seam. Only need y's
    public int[] findHorizontalSeam() {               
        // only maintain one line is enough   
        distTo = new double[height()];
        // pixTo[x][y] return the y only, because x must be x - 1
        pixTo = new short[width()][height()];
        // initialize first row, same as connecting them to top virtual node, then relax them.
        // don't need the pixTo[i][0] information.
        for (int j = 0; j < height(); j++) {
            distTo[j] = energy[0][j];
        }
        
            
        // visit vertices in toplogical order, use virtual node at left
        // do nothing for last column           
        for (int i = 0; i < width() - 1; i++) {                
            double distTemp, distNext;
            // for j = 0, initialize right and down-right, relax all
            distTemp = distTo[0];
            distTo[0] = Double.POSITIVE_INFINITY;
            relaxH(i + 1, 0, 0, distTemp);                
            // store current distNext
            distNext = distTo[1];
            distTo[1] = Double.POSITIVE_INFINITY;
            relaxH(i + 1, 1, 0, distTemp); 
                
            // for 0 < j < height() - 1, initialize the down-right only, relax all
            for (int j = 1; j < height() - 1; j++) { 
                distTemp = distNext;
                relaxH(i + 1, j, j, distTemp);
                relaxH(i + 1, j - 1, j, distTemp);
                   
                distNext = distTo[j + 1];
                distTo[j + 1] = Double.POSITIVE_INFINITY;
                relaxH(i + 1, j + 1, j, distTemp);                     
            }
            // for j = height() - 1, don't initialize, just relax all 
            distTemp = distNext;
            relaxH(i + 1, height() - 1, height() - 1, distTemp);
            relaxH(i + 1, height() - 2, height() - 1, distTemp);
        }
        // done, find the minSeam
        double minDist = Double.POSITIVE_INFINITY;        
        int minPix = 0;
        for  (int j = 0; j < height() - 1; j++) {
            if (distTo[j] < minDist) {                    
                minDist = distTo[j];
                minPix = j;
            }            
        }
        // store minSeam    
        int[] minSeam = new int[width()];
        int from = minPix;
        minSeam[width() - 1] = from;
        for (int i = width() - 1; i > 0; i--) {
            // store the y only   
            from = pixTo[i][from];         
            minSeam[i - 1] = from;             
        }                
        //StdOut.printf("%d %d\n", toV(0, y), from);
        return minSeam;
    }

    // sequence of indices for vertical seam. Only need x's 
    public int[] findVerticalSeam() {
        // only maintain one line is enough     
        distTo = new double[width()];
        // pixTo returns the x only, because y must be y - 1
        pixTo = new short[width()][height()];
        
        // initialize first row, same as connecting them to top virtual node, then relax them.
        // don't need the pixTo[i][0] information.
        for (int i = 0; i < width(); i++) {
            distTo[i] = energy[i][0];
        }
            
        // visit vertices in toplogical order, use virtual node at top
        // do nothing for last row
        for (int j = 0; j < height() - 1; j++) {              
            double distTemp, distNext;
            // for i = 0, initialize below and down-right, relax all
            distTemp = distTo[0];
            distTo[0] = Double.POSITIVE_INFINITY;
            relaxV(0, j + 1, 0, distTemp);

            distNext = distTo[1];
            distTo[1] = Double.POSITIVE_INFINITY;
            relaxV(1, j + 1, 0, distTemp); 
            // for 0 < i < width() - 1, initialize the down-right only, relax all
            for (int i = 1; i < width() - 1; i++) {                
                distTemp = distNext;
                relaxV(i, j + 1, i, distTemp);
                relaxV(i - 1, j + 1, i, distTemp);
               
                distNext = distTo[i + 1];
                distTo[i + 1] = Double.POSITIVE_INFINITY;
                relaxV(i + 1, j + 1, i, distTemp); 
            }
            // for i = width() - 1, don't initialize, just relax all 
            distTemp = distNext;
            relaxV(width() - 1, j + 1, width() - 1, distTemp);
            relaxV(width() - 2, j + 1, width() - 1, distTemp);
            
        }
    
        // done, find the minSeam  
        double minDist = Double.POSITIVE_INFINITY;
        int minPix = 0;
        for  (int i = 0; i < width(); i++) {
            if (distTo[i] < minDist) {                    
                minDist = distTo[i];
                minPix = i;
            }
        }           
        
        // store minSeam           
        int[] minSeam = new int[height()];            
        int from = minPix;
        minSeam[height() - 1] = from; 
        for (int j = height() - 1; j > 0; j--) {
            // store the x only            
            from = pixTo[from][j];
            minSeam[j - 1] = from; 
        }                
        //StdOut.println(minSeam[1]);
        return minSeam;
    } 

    // remove horizontal seam from current picture
    public void removeHorizontalSeam(int[] seam) {
        if (seam == null) throw new NullPointerException("No seam!");
        if (seam.length != width()) throw new IllegalArgumentException("Not a valid horizontal seam!");  
        if (height() < 2) throw new IllegalArgumentException("Can't remove!");
        hc -= 1;
        for (int i = 0; i < width(); i++) {
            if (i != 0 && (seam[i] > seam[i - 1] + 1 || seam[i] < seam[i - 1] -1)) throw new IllegalArgumentException("Not a valid seam!");
            if (seam[i] < 0 || seam[i] > height()) throw new IllegalArgumentException("Not a valid seam!");
            // shift all the way up
            
            for (int j = seam[i]; j < height(); j++) {
                temp[j][i] = temp[j + 1][i];
                //ex(temp, j, i, j + 1, i);               
            }                        
        }
        
        // change the energy matrix. Skip first and last colunm
        for (int i = 1; i < width() - 1; i++) {
            if (seam[i] == height()) { 
                energy[i][seam[i] - 1] = 195075; 
                continue;
            }
            energy[i][seam[i]] = (int) energy(i, seam[i]);
            if (seam[i] != 0) energy[i][seam[i] - 1] = (int) energy(i, seam[i] - 1);
            if (seam[i] + 1 < height()) {
                /*
                for (int j = seam[i] + 1; j < height(); j++) {
                    energy[i][j] = energy[i][j + 1];             
                }
                */
                System.arraycopy(energy[i], seam[i] + 2, energy[i], seam[i] + 1, height() - seam[i] - 1);
            }
            
        }
        
        return;
    }

    // remove vertical seam from current picture 
    public void removeVerticalSeam(int[] seam) {
        if (seam == null) throw new NullPointerException("No seam!");
        if (seam.length != height()) throw new IllegalArgumentException("Not a valid vertical seam!");   
        if (width() < 2) throw new IllegalArgumentException("Can't remove!");
        wc -= 1;
        for (int j = 0; j < height(); j++) {
            if (j != 0 && (seam[j] > seam[j - 1] + 1 || seam[j] < seam[j - 1] -1)) throw new IllegalArgumentException("Not a valid seam!");
            if (seam[j] < 0 || seam[j] > width()) throw new IllegalArgumentException("Not a valid seam!");
            // shift all the way left
            /*
            for (int i = seam[j]; i < width(); i++) {
                temp[j][i] = temp[j][i + 1];            
            }
            */
            System.arraycopy(temp[j], seam[j] + 1, temp[j], seam[j], width() - seam[j]);
        }
        
        // change the energy matrix. Skip first and last row
        for (int j = 1; j < height() - 1; j++) {
            if (seam[j] == width()) { 
                energy[seam[j] - 1][j] = 195075; 
                continue;
            }
            energy[seam[j]][j] = (int) energy(seam[j], j);
            if (seam[j] != 0) energy[seam[j] - 1][j] = (int) energy(seam[j] - 1, j);
            if (seam[j] + 1 < width()) {
                for (int i = seam[j] + 1; i < width(); i++) {
                    energy[i][j] = energy[i + 1][j];
                    //exe(energy, toV(i, j), toV(i + 1, j));               
                }         
            }
        }
        
        return;
    } 

    // relax H or V
    // xfrom = x - 1
    private void relaxH(int x, int y, int yfrom, double dist) {
        if (distTo[y] > dist + energy[x][y]) {
            distTo[y] = dist + energy[x][y];
            pixTo[x][y] = (short) yfrom;
        } 
    }
    // yfrom = y - 1
    private void relaxV(int x, int y, int xfrom, double dist) {       
        if (distTo[x] > dist + energy[x][y]) {
            distTo[x] = dist + energy[x][y];
            pixTo[x][y] = (short) xfrom;
        } 
    }
    // convert ARGB number into RGB
    private int toR(int argb) { return (argb >> 16) & 0xFF; } 
    private int toG(int argb) { return (argb >> 8) & 0xFF; } 
    private int toB(int argb) { return (argb) & 0xFF; }    

    // test client
    public static void main(String[] args) { 
        /*
        Picture inputImg = new Picture(args[0]);
        System.out.printf("image is %d pixels wide by %d pixels high.\n", inputImg.width(), inputImg.height());
        
        SeamCarver sc = new SeamCarver(inputImg);
        
        
        PrintSeams ps = new PrintSeams();
        ps.printHorizontalSeam(sc);
        ps.printVerticalSeam(sc);        
       
        
        System.out.printf("Printing energy calculated for each pixel.\n");      

        for (int j = 0; j < sc.height(); j++)
        {
            for (int i = 0; i < sc.width(); i++)
                System.out.printf("%9.0f ", sc.energy(i, j));
            System.out.println();
        }

        int[] vs = sc.findVerticalSeam();
        sc.removeVerticalSeam(vs);
        
        System.out.println("Energy after removing one vertical seam:");      
        for (int j = 0; j < sc.height(); j++)
        {
            for (int i = 0; i < sc.width(); i++)
                System.out.printf("%9.0f ", sc.energy(i, j));
            System.out.println();
        }
        
        //System.out.println("After removing one vertical seam, find seam again:");      
        //ps.printVerticalSeam(sc);
        //sc.energy(12, 11);
        int[] right = new int[sc.height()];
        for (int i = 0; i < sc.height(); i++) {
            right[i] = sc.width() - 1;
        }
        sc.removeVerticalSeam(right);

        System.out.println("Energy after removing right-most vertical seam:");      
        for (int j = 0; j < sc.height(); j++)
        {
            for (int i = 0; i < sc.width(); i++)
                System.out.printf("%9.0f ", sc.energy(i, j));
            System.out.println();
        }

        vs = sc.findVerticalSeam();
        sc.removeVerticalSeam(vs);

        System.out.println("Energy after removing one vertical seam:");      
        for (int j = 0; j < sc.height(); j++)
        {
            for (int i = 0; i < sc.width(); i++)
                System.out.printf("%9.0f ", sc.energy(i, j));
            System.out.println();
        }
        */
    }
}
