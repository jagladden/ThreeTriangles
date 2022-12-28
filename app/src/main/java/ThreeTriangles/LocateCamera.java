package ThreeTriangles;

import java.util.ArrayList;
import java.util.HashMap;
import static java.lang.Math.*;

public class LocateCamera {
    ArrayList<DetectedTag> tags = new ArrayList<DetectedTag>();
    HashMap<Integer, LandMark> landmarks = new HashMap<Integer, LandMark>();
    ArrayList<TrianglePair> trianglePairs = new ArrayList<TrianglePair>();
    ArrayList<Result> results = new ArrayList<Result>();

    public ArrayList<TrianglePair> solve() {
        
        //Sort by increasing angle
        tags.sort(null);
        int nTags = tags.size();

        //Add the associate landmark to each tag
        for (int i = 0; i < nTags; ++i) {
            tags.get(i).lm = landmarks.get(tags.get(i).id); 
        }

        //Build the list of possible triangle pairs with adjoining interior sides
        //formed by the observer positions and the observed tags
        for (int i = 0; i < nTags; ++i) {
            for (int j = i+1; j < nTags; ++j) {
                for(int k = j+1; k < nTags; ++k) {
                    TrianglePair trp = new TrianglePair(tags.get(i), tags.get(j), tags.get(k));
                    trp.computeSensitivity(.01);
                    trianglePairs.add(trp); 
                }
              
            }
        }

        //The constructor for the TrianglePair object already computed the observer position for each
        //set of of three observed landmarks. 

        return trianglePairs;
    }

    public void addDetectedTag(DetectedTag dt) {
        tags.add(dt);
    }

    public void clearTags() {
        tags.clear();
    }

    public void addLandmark(LandMark lm) {
        landmarks.put(lm.tagID, lm);
    }

    public static class LandMark {
        public Point p;
        public int tagID;

        public LandMark(Point p, int tagID) {
            this.p = p;
            this.tagID = tagID;
        }

    }

    public static class DetectedTag implements Comparable<DetectedTag> {
        int id;
        public double angle;
        public LandMark lm;
        

        public DetectedTag (int id, double angle) {
            this.id = id;
            this.angle = angle;
        }

        //Sort by angle
        @Override
        public int compareTo(DetectedTag o) {
            if (angle < o.angle) {
                return -1;
            } else if (angle > o.angle) {
                return 1;
            } else {
                return 0;
            }
        }
    }

    public static class Point {
        public double x;
        public double y;

        public Point(double x, double y) {
            this.x = x;
            this.y = y;
        }

        double distance(Point p) {
            return sqrt((p.x-this.x) * (p.x-this.x) +  (p.y-this.y) * (p.y-this.y));
        }
    }

    public static class Result {
        TrianglePair trp;
        Point p;
         
        public Result(TrianglePair trp, Point p) {
            this.trp = trp;
            this.p = p;
        }
    }

    public static class TrianglePair {
        DetectedTag tag1;
        DetectedTag tag2;
        DetectedTag tag3;
        double a;
        double b;
        double e;
        double f;
        double C1;
        double C2;
        double C3;
        double rotationFromWCS;
        Point observer;
        double errorSensitivity;
      
        public TrianglePair(DetectedTag tag1, DetectedTag tag2, DetectedTag tag3) {
            this.tag1 = tag1;
            this.tag2 = tag2;
            this.tag3 = tag3;

            C1 = tag2.lm.p.distance(tag1.lm.p);
            C2 = tag3.lm.p.distance(tag2.lm.p);
            C3 = tag3.lm.p.distance(tag1.lm.p);

            //Compute the rotation angle that would be required to place F1 on theY-axis
            //of our local coordinate system assuming that translation had already occurred
            //such that F3 is at the origin.
            rotationFromWCS = PI/2 - atan2(tag1.lm.p.y - tag3.lm.p.y, tag1.lm.p.x - tag3.lm.p.x);

            //Now rotate F2 to determine it's X coordinate. If the result is positive, then angles
            //e and f are considered to be positive for purposes of this algorithm, otherwise
            double f2x = (tag2.lm.p.x - tag3.lm.p.x) * cos(rotationFromWCS)
                       - (tag2.lm.p.y - tag3.lm.p.y) * sin(rotationFromWCS);

            e = acos((square(C2) + square(C3) - square(C1)) / (2 * C2 * C3));
            f = acos((square(C1) + square(C3) - square(C2)) / (2 * C1 * C3));
            if (f2x < 0) {
                e = -e;
                f = -f;
            }

            observer = computeObserverPosition(0, 0);
        }

        public Point computeObserverPosition(double ha, double hb) {
            a = tag2.angle - tag1.angle + ha;
            b = tag3.angle - tag2.angle + hb;
           
            double cRatio1 = C1 / sin(a);
            double cRatio2 = C2 / sin(b);
            double k = a + b + e + f;

            //Compute interior angle c (opposite side A) of triangle A-B2-C2
            double num = cRatio1 * sin(k);
            double denom = cRatio1*cos(k) - cRatio2;
            double z = atan2(num, denom);

            double c;
            if (z > 0) {
                c = PI - z;
            } else {
                c = -z;
            }

            //Now, from c, get g and h
            double g = PI - (b + c);
            double h = PI/2 - (e + c);

            //And hypotenuse B2
            double B2 = cRatio2 * sin(g);

            //The observer position in the local coordinate system
            Point pLocal = new Point(B2 * cos(h), B2 * sin(h));
      
            //The model for which the observer X,Y were calculated uses a coordinate system
            //in which tag1 and tag3 both lie on the X axis, and tag3 is at the origin. Thus
            //the tags were rotated by an amount b = pi/2 - a, where a is the angle of the line
            //from tag3 to tag 1 relative to the WCS X axis.  Thus, the computed X, Y postion
            //needs to be rotated by angle -b before translation.
            double ra = -rotationFromWCS;
            double xr = pLocal.x * cos(ra) - pLocal.y * sin(ra);
            double yr = pLocal.y * cos(ra) + pLocal.x * sin(ra);
            double xrt = xr + tag3.lm.p.x;
            double yrt = yr + + tag3.lm.p.y;
            
            return new Point(xrt, yrt);
        }

        //Use a four point derivative formula to estimate the sensitivity of the computed observer
        //position to a small perbutations in measured angles 'a' and 'b'. For the purposes of the
        //derivative the 'position' is the absolute distance from the origin. Sum the two
        //derivatives to create our sensitivity metric.
        public void computeSensitivity(double delta) {
            Point p0 = new Point(0, 0);
            double[] coef = {-11, 18, -9, 2};
            double mult = 1/(6 * delta);

            double sum_a = 0;
            double sum_b = 0;
            for (int i = 0; i < 4; ++i) {
                Point p1 = computeObserverPosition(delta * i, 0);
                sum_a += coef[i] * p1.distance(p0);
                Point p2 = computeObserverPosition(0, delta * i);
                sum_b += coef[i] * p2.distance(p0);
            }
            errorSensitivity = (sum_a + sum_b) * mult;
        }

    }

    public static double square(double x) {
        return x*x;
    }
    
}
