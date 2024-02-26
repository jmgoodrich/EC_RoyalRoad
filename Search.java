/******************************************************************************
*  A Teaching GA					  Developed by Hal Stringer & Annie Wu, UCF
*  Version 2, January 18, 2004
*******************************************************************************/

import java.io.*;
import java.util.*;
import java.text.*;

public class Search {

/*******************************************************************************
*                           INSTANCE VARIABLES                                 *
*******************************************************************************/

/*******************************************************************************
*                           STATIC VARIABLES                                   *
*******************************************************************************/

	public static FitnessFunction problem;

	public static Chromo[] member;
	public static Chromo[] child;

	public static Chromo bestOfGenChromo;
	public static int bestOfGenR;
	public static int bestOfGenG;
	public static Chromo bestOfRunChromo;
	public static int bestOfRunR;
	public static int bestOfRunG;
	public static Chromo bestOverAllChromo;
	public static int bestOverAllR;
	public static int bestOverAllG;

	public static double sumRawFitness;
	public static double sumRawFitness2;	// sum of squares of fitness
	public static double sumSclFitness;
	public static double sumProFitness;
	public static double defaultBest;
	public static double defaultWorst;

	public static double averageRawFitness;
	public static double stdevRawFitness;

	public static int G;
	public static int R;
	public static Random r = new Random();
	private static double randnum;

	private static int memberIndex[];
	private static double memberFitness[];
	private static int TmemberIndex;
	private static double TmemberFitness;

	private static double fitnessStats[][];  // 0=Avg, 1=Best

	public static int numBlocks[][];
	public static double avgNumBlocks[][];
	private static double stdError;
	public static int extraBlocks;

	public static double avgNumGenerations;
	public static double medianNumGenerations;
	public static int[] numGenerations;
	public static int worstGeneration;

	private static boolean foundSolution;
	private static int solution;

/*******************************************************************************
*                              CONSTRUCTORS                                    *
*******************************************************************************/


/*******************************************************************************
*                             MEMBER METHODS                                   *
*******************************************************************************/


/*******************************************************************************
*                             STATIC METHODS                                   *
*******************************************************************************/

	public static String formatDouble(double value) {
		if (value == (long) value) {
			return String.format("%d", (long) value);
		} else {
			return String.format("%.4f", value).replaceAll("0*$", "").replaceAll("\\.$", "");
		}
	}

	public static void main(String[] args) throws java.io.IOException{

		Calendar dateAndTime = Calendar.getInstance(); 
		Date startTime = dateAndTime.getTime();

	//  Read Parameter File
		System.out.println("\nParameter File Name is: " + args[0] + "\n");
		Parameters parmValues = new Parameters(args[0]);

	//  Write Parameters To Summary Output File
		String summaryFileName = Parameters.expID + "_summary.txt";
		String blockFileName = Parameters.expID + "_blocks_summary.txt";
		FileWriter summaryOutput = new FileWriter(summaryFileName);
		FileWriter blockSummaryOutput = new FileWriter(blockFileName);
		parmValues.outputParameters(summaryOutput);

	//	Set up Fitness Statistics matrix
		fitnessStats = new double[3][Parameters.generations];
		for (int i=0; i<Parameters.generations; i++){
			fitnessStats[0][i] = 0;
			fitnessStats[1][i] = 0;
			fitnessStats[2][i] = 0;
		}

	//	Set up building blocks statistics
		extraBlocks = 0;
		if(Parameters.isNonlinear) {
			int j = Parameters.numGenes * 2;
			while(j < Parameters.geneSize * Parameters.numGenes) {
				extraBlocks += Parameters.geneSize * Parameters.numGenes / j;
				j *= 2;
			}
		}
		avgNumBlocks = new double[Parameters.generations][Parameters.numGenes + extraBlocks];	//	Sets all values to 0
		avgNumGenerations = 0;
		numGenerations = new int[Parameters.numRuns];
		worstGeneration = 0;

	//	Problem Specific Setup - For new new fitness function problems, create
	//	the appropriate class file (extending FitnessFunction.java) and add
	//	an else_if block below to instantiate the problem.
 
		if (Parameters.problemType.equals("NM")){
				problem = new NumberMatch();
		}
		else if (Parameters.problemType.equals("RR")){
			problem = new RoyalRoad();
		}
		else System.out.println("Invalid Problem Type");

		System.out.println(problem.name);

	//	Initialize RNG, array sizes and other objects
		r.setSeed(Parameters.seed);
		memberIndex = new int[Parameters.popSize];
		memberFitness = new double[Parameters.popSize];
		member = new Chromo[Parameters.popSize];
		child = new Chromo[Parameters.popSize];
		bestOfGenChromo = new Chromo();
		bestOfRunChromo = new Chromo();
		bestOverAllChromo = new Chromo();
		if(Parameters.isNonlinear)
			solution = Parameters.geneSize * Parameters.numGenes * (int)(Math.log(Parameters.numGenes) / Math.log(2));
		else 
			solution = Parameters.geneSize * Parameters.numGenes;

		if (Parameters.minORmax.equals("max")){
			defaultBest = 0;
			defaultWorst = 999999999999999999999.0;
		}
		else{
			defaultBest = 999999999999999999999.0;
			defaultWorst = 0;
		}

		bestOverAllChromo.rawFitness = defaultBest;

		//  Start program for multiple runs
		for (R = 1; R <= Parameters.numRuns; R++){

			bestOfRunChromo.rawFitness = defaultBest;
			System.out.println();

			//	Initialize First Generation
			for (int i=0; i<Parameters.popSize; i++){
				member[i] = new Chromo();
				child[i] = new Chromo();
			}

			//	Initialize schema info
			numBlocks = new int[Parameters.generations][Parameters.numGenes + extraBlocks];
			foundSolution = false;

			//	Begin Each Run
			for (G=0; G<Parameters.generations; G++){

				sumProFitness = 0;
				sumSclFitness = 0;
				sumRawFitness = 0;
				sumRawFitness2 = 0;
				bestOfGenChromo.rawFitness = defaultBest;

				// Schema information
				numBlocks[G] = new int[Parameters.numGenes + extraBlocks];

				//	Test Fitness of Each Member
				for (int i=0; i<Parameters.popSize; i++){

					member[i].rawFitness = 0;
					member[i].sclFitness = 0;
					member[i].proFitness = 0;

					problem.doRawFitness(member[i]);

					//	Update building block statistics
					for(int j = 0; j < member[i].hasBlock.length; j++) 
						if(member[i].hasBlock[j]) numBlocks[G][j]++;

					sumRawFitness = sumRawFitness + member[i].rawFitness;
					sumRawFitness2 = sumRawFitness2 +
						member[i].rawFitness * member[i].rawFitness;

					if (Parameters.minORmax.equals("max")){
						if (member[i].rawFitness > bestOfGenChromo.rawFitness){
							Chromo.copyB2A(bestOfGenChromo, member[i]);
							bestOfGenR = R;
							bestOfGenG = G;
						}
						if (member[i].rawFitness > bestOfRunChromo.rawFitness){
							Chromo.copyB2A(bestOfRunChromo, member[i]);
							bestOfRunR = R;
							bestOfRunG = G;
						}
						if (member[i].rawFitness > bestOverAllChromo.rawFitness){
							Chromo.copyB2A(bestOverAllChromo, member[i]);
							bestOverAllR = R;
							bestOverAllG = G;
						}
					}
					else {
						if (member[i].rawFitness < bestOfGenChromo.rawFitness){
							Chromo.copyB2A(bestOfGenChromo, member[i]);
							bestOfGenR = R;
							bestOfGenG = G;
						}
						if (member[i].rawFitness < bestOfRunChromo.rawFitness){
							Chromo.copyB2A(bestOfRunChromo, member[i]);
							bestOfRunR = R;
							bestOfRunG = G;
						}
						if (member[i].rawFitness < bestOverAllChromo.rawFitness){
							Chromo.copyB2A(bestOverAllChromo, member[i]);
							bestOverAllR = R;
							bestOverAllG = G;
						}
					}
				}

				// Accumulate fitness statistics
				fitnessStats[0][G] += sumRawFitness / Parameters.popSize;
				fitnessStats[1][G] += bestOfGenChromo.rawFitness;

				// Found a solution?
				if(bestOfGenChromo.rawFitness == solution) foundSolution = true;

				averageRawFitness = sumRawFitness / Parameters.popSize;
				stdevRawFitness = Math.sqrt(
							Math.abs(sumRawFitness2 - 
							sumRawFitness*sumRawFitness/Parameters.popSize)
							/
							(Parameters.popSize-1)
							);
				fitnessStats[2][G] += stdevRawFitness;

				// Output generation statistics to screen
				System.out.println(R + "\t" + G +  "\t" + (int)bestOfGenChromo.rawFitness + "\t" + averageRawFitness + "\t" + stdevRawFitness);

				// Output generation statistics to summary file
				summaryOutput.write(" R ");
				Hwrite.right(R, 3, summaryOutput);
				summaryOutput.write(" G ");
				Hwrite.right(G, 3, summaryOutput);
				Hwrite.right((int)bestOfGenChromo.rawFitness, 7, summaryOutput);
				Hwrite.right(averageRawFitness, 11, 3, summaryOutput);
				Hwrite.right(stdevRawFitness, 11, 3, summaryOutput);
				summaryOutput.write("\n");

				// Output generation statistics to block file
				//blockSummaryOutput.write(formattedR + "" + formattedG);
				//for(int i = 0; i < numBlocks.length; i++) {
				//	String formatted = String.format("%9d", numBlocks[G][i]);
				//	blockSummaryOutput.write(formatted + " ");
				//}
				//blockSummaryOutput.write("\n");


		// *********************************************************************
		// **************** SCALE FITNESS OF EACH MEMBER AND SUM ***************
		// *********************************************************************

				switch(Parameters.scaleType){

				case 0:     // No change to raw fitness
					for (int i=0; i<Parameters.popSize; i++){
						member[i].sclFitness = member[i].rawFitness + .000001;
						sumSclFitness += member[i].sclFitness;
					}
					break;

				case 1:     // Fitness not scaled.  Only inverted.
					for (int i=0; i<Parameters.popSize; i++){
						member[i].sclFitness = 1/(member[i].rawFitness + .000001);
						sumSclFitness += member[i].sclFitness;
					}
					break;

				case 2:     // Fitness scaled by Rank (Maximizing fitness)

					//  Copy genetic data to temp array
					for (int i=0; i<Parameters.popSize; i++){
						memberIndex[i] = i;
						memberFitness[i] = member[i].rawFitness;
					}
					//  Bubble Sort the array by floating point number
					for (int i=Parameters.popSize-1; i>0; i--){
						for (int j=0; j<i; j++){
							if (memberFitness[j] > memberFitness[j+1]){
								TmemberIndex = memberIndex[j];
								TmemberFitness = memberFitness[j];
								memberIndex[j] = memberIndex[j+1];
								memberFitness[j] = memberFitness[j+1];
								memberIndex[j+1] = TmemberIndex;
								memberFitness[j+1] = TmemberFitness;
							}
						}
					}
					//  Copy ordered array to scale fitness fields
					for (int i=0; i<Parameters.popSize; i++){
						member[memberIndex[i]].sclFitness = i;
						sumSclFitness += member[memberIndex[i]].sclFitness;
					}

					break;

				case 3:     // Fitness scaled by Rank (minimizing fitness)

					//  Copy genetic data to temp array
					for (int i=0; i<Parameters.popSize; i++){
						memberIndex[i] = i;
						memberFitness[i] = member[i].rawFitness;
					}
					//  Bubble Sort the array by floating point number
					for (int i=1; i<Parameters.popSize; i++){
						for (int j=(Parameters.popSize - 1); j>=i; j--){
							if (memberFitness[j-i] < memberFitness[j]){
								TmemberIndex = memberIndex[j-1];
								TmemberFitness = memberFitness[j-1];
								memberIndex[j-1] = memberIndex[j];
								memberFitness[j-1] = memberFitness[j];
								memberIndex[j] = TmemberIndex;
								memberFitness[j] = TmemberFitness;
							}
						}
					}
					//  Copy array order to scale fitness fields
					for (int i=0; i<Parameters.popSize; i++){
						member[memberIndex[i]].sclFitness = i;
						sumSclFitness += member[memberIndex[i]].sclFitness;
					}

					break;

				case 4:		//	Fitness scaled using sigma scaling

					for (int i=0; i<Parameters.popSize; i++){

						// SD = 0: All same scl fitness
						if(stdevRawFitness == 0.0) {
							member[i].sclFitness = 1;
						}

						// Sigma scaling
						else {
							member[i].sclFitness = 1 + ((member[i].rawFitness - averageRawFitness) / (2.0 * stdevRawFitness));
							if(member[i].sclFitness > 1.5) member[i].sclFitness = 1.5;
						}
						sumSclFitness += member[i].sclFitness;
					}
					
					break;

				default:
					System.out.println("ERROR - No scaling method selected");
				}


		// *********************************************************************
		// ****** PROPORTIONALIZE SCALED FITNESS FOR EACH MEMBER AND SUM *******
		// *********************************************************************

				for (int i=0; i<Parameters.popSize; i++){
					member[i].proFitness = member[i].sclFitness/sumSclFitness;
					sumProFitness = sumProFitness + member[i].proFitness;
				}

		// *********************************************************************
		// ************ CROSSOVER AND CREATE NEXT GENERATION *******************
		// *********************************************************************

				int parent1 = -1;
				int parent2 = -1;

				//  Assumes always two offspring per mating
				for (int i=0; i<Parameters.popSize; i=i+2){

					//	Select Two Parents
					parent1 = Chromo.selectParent();
					parent2 = parent1;
					while (parent2 == parent1){
						parent2 = Chromo.selectParent();
					}

					//	Crossover Two Parents to Create Two Children
					randnum = r.nextDouble();
					if (randnum < Parameters.xoverRate){
						Chromo.mateParents(parent1, parent2, member[parent1], member[parent2], child[i], child[i+1]);
					}
					else {
						Chromo.mateParents(parent1, member[parent1], child[i]);
						Chromo.mateParents(parent2, member[parent2], child[i+1]);
					}
				} // End Crossover

				//	Mutate Children
				for (int i=0; i<Parameters.popSize; i++){
					child[i].doMutation();
				}

				//	Swap Children with Last Generation
				for (int i=0; i<Parameters.popSize; i++){
					Chromo.copyB2A(member[i], child[i]);
				}

				// Exit condition
				if(foundSolution) break;

			} //  Repeat the above loop for each generation or until a solution is found

			Hwrite.left(bestOfRunR, 4, summaryOutput);
			Hwrite.right(bestOfRunG, 4, summaryOutput);

			problem.doPrintGenes(bestOfRunChromo, summaryOutput);

			System.out.println(R + "\t" + "B" + "\t"+ (int)bestOfRunChromo.rawFitness);

			// Schema statistics
			if(G == Parameters.generations) G--;
			numGenerations[R - 1] = G;
			for(int i = 0; i <= G; i++) {
				for(int j = 0; j < Parameters.numGenes + extraBlocks; j++) {
					avgNumBlocks[i][j] += numBlocks[i][j];
				}
			}
			
			if(worstGeneration < G) worstGeneration = G;

		} //End of a Run

		// Update schema information
		for(int i = 0; i < Parameters.numRuns; i++) avgNumGenerations += numGenerations[i];
		avgNumGenerations /= 1.0 * Parameters.numRuns;

		Arrays.sort(numGenerations);
		medianNumGenerations = Parameters.numRuns % 2 == 0 ? 
			(numGenerations[Parameters.numRuns / 2 - 1] + numGenerations[Parameters.numRuns / 2]) / 2.0 : 
			numGenerations[Parameters.numRuns / 2];

		for(int i = 0; i <= worstGeneration; i++) {
			for(int j = 0; j < Parameters.numGenes + extraBlocks; j++) {
				avgNumBlocks[i][j] /= (1.0 * Parameters.numRuns);
			}
		}

		stdError = 0.0;
		for(int i = 0; i < Parameters.numRuns; i++)
			stdError += Math.pow(numGenerations[i] - avgNumGenerations, 2);
		stdError /= ((1.0 * Parameters.numRuns) - 1);
		stdError = Math.sqrt(stdError);
		stdError /= Math.sqrt(Parameters.numRuns);


		Hwrite.left("B", 8, summaryOutput);

		problem.doPrintGenes(bestOverAllChromo, summaryOutput);

		//	Output Fitness Statistics matrix
		summaryOutput.write("Gen            AvgFit              BestFit             StdDev\n");
		for (int i=0; i<=worstGeneration; i++){
			Hwrite.left(i, 15, summaryOutput);
			Hwrite.left(fitnessStats[0][i]/Parameters.numRuns, 20, 2, summaryOutput);
			Hwrite.left(fitnessStats[1][i]/Parameters.numRuns, 20, 2, summaryOutput);
			Hwrite.left(fitnessStats[2][i]/Parameters.numRuns, 20, 2, summaryOutput);
			summaryOutput.write("\n");
		}

		summaryOutput.write("\n");
		summaryOutput.close();


		//	Output block statistics
		Parameters.outputParameters(blockSummaryOutput);
		blockSummaryOutput.write("\nNumber of schemas: " + (Parameters.numGenes + extraBlocks) + "\n---------------------\n");
		for(int i = 1; i <= Parameters.numGenes; i++) blockSummaryOutput.write("s" + i + ": " + Parameters.geneSize + "\n");
		if(Parameters.isNonlinear) {
			int k = Parameters.numGenes + 1;
			int div = 2;
			while(div < Parameters.numGenes) {
				// Current fitness weight
				int val = Parameters.geneSize * div;
				for(int i = 0; i < Parameters.numGenes / div; i++)
					blockSummaryOutput.write("s" + (k++) + ": " + val + "\n");
				div *= 2;
			}
		}

		blockSummaryOutput.write("\n\nAverage # of Gens = " + avgNumGenerations);
		blockSummaryOutput.write("\nAverage # of Gens Standard Error = " + stdError);
		blockSummaryOutput.write("\nMedian # of Gens = " + medianNumGenerations);

		blockSummaryOutput.write("\n\nAverage # of Schema Per Generation");
		blockSummaryOutput.write("\nGen    ");
		for(int i = 1; i <= (Parameters.numGenes + extraBlocks); i++) {
			String formatted = String.format("s" + "%-8d", i);
			blockSummaryOutput.write(formatted);
		}
		blockSummaryOutput.write("\n");

		for(int i = 0; i <= worstGeneration; i++) {
			String formatted = String.format("%-7d", i);
			blockSummaryOutput.write(formatted);
			for(int j = 0; j < avgNumBlocks[i].length; j++) {
				formatted = String.format("%-8s", formatDouble(avgNumBlocks[i][j]));
				blockSummaryOutput.write(formatted + " ");
			}
			blockSummaryOutput.write("\n");
		}

		blockSummaryOutput.write("\n");
		blockSummaryOutput.close();

		System.out.println();
		System.out.println("Start:  " + startTime);
		dateAndTime = Calendar.getInstance(); 
		Date endTime = dateAndTime.getTime();
		System.out.println("End  :  " + endTime);

		System.out.println("Solution = " + solution);
		System.out.println("Worst gen = " + worstGeneration);

	} // End of Main Class

}   // End of Search.Java ******************************************************

