import java.io.*;
import java.util.*;
import java.text.*;

public class OneMax extends FitnessFunction{

/*******************************************************************************
*                            INSTANCE VARIABLES                                *
*******************************************************************************/


/*******************************************************************************
*                            STATIC VARIABLES                                  *
*******************************************************************************/


/*******************************************************************************
*                              CONSTRUCTORS                                    *
*******************************************************************************/

	public RoyalRoad(){
		name = "RoyalRoad Problem";
	}

/*******************************************************************************
*                                MEMBER METHODS                                *
*******************************************************************************/

//  COMPUTE A CHROMOSOME'S RAW FITNESS *************************************

	public void doRawFitness(Chromo X){

		X.rawFitness = 0;
        int steps = Parameters.intermediateBlocks;


		for(int i = 0; i < Parameters.numGenes; i++) {
			boolean flag = false;
			for(int j = 0; j < Parameters.geneSize; j++) {

				// Not a block
				if(X.chromo.charAt(i * Parameters.geneSize + j) == '0') {
					flag = true;
					break;
				}
			}

			// Is this a block?
			if(!flag) {
				X.hasBlock[i] = true;
				X.rawFitness += Parameters.geneSize;
			}
		}

        // if(steps > 0){
        //     for(int k = steps; k > 0; k--){
        //         for(int m = 0; m < Parameters.numGenes; m = m+Math.pow(2, k)){
        //             for(int n = 0; n < Parameters.numGenes/steps; n++) {
        //                 if(!X.hasBlock[n+m]){
        //                     break;
        //                 }
        //                 X.rawFitness += Parameters.geneSize * Math.pow(2, k);
        //             }
        //         }
        //     }
        // }
	}

//  PRINT OUT AN INDIVIDUAL GENE TO THE SUMMARY FILE *********************************

	public void doPrintGenes(Chromo X, FileWriter output) throws java.io.IOException{

		for (int i=0; i<Parameters.numGenes; i++){
			Hwrite.right(X.getGeneAlpha(i),11,output);
		}
		output.write("   RawFitness");
		output.write("\n        ");
		for (int i=0; i<Parameters.numGenes; i++){
			Hwrite.right(X.getPosIntGeneValue(i),11,output);
		}
		Hwrite.right((int) X.rawFitness,13,output);
		output.write("\n\n");
		return;
	}

/*******************************************************************************
*                             STATIC METHODS                                   *
*******************************************************************************/

}   // End of RoyalRoad.java ******************************************************