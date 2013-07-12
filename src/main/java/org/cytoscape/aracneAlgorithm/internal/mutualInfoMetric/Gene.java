package org.cytoscape.aracneAlgorithm.internal.mutualInfoMetric;


public  class Gene {
        public double x;
        public int xi; // index of x
        public int maId;
        
        public Gene(int mId, double value) {
            maId = mId;
            x = value;
        }
        
        @Override
        public boolean equals(Object o) {
            Gene gp = (Gene)o;
            return (gp.maId == maId);
        }

        @Override
        public int hashCode() {
            int hash = 7;
            hash = 23 * hash + this.maId;
            return hash;
        }
        
        @Override
        public String toString() {
            return "x:" + x + "xi:" + x + "maId:" + maId;
        }
}
