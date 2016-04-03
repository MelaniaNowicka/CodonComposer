using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Linq;
using System.Text;
using System.Threading;
using System.Threading.Tasks;

namespace CodonOptimizer.Classes
{
    /// <summary>
    /// A class for optimization of ORF
    /// includes implementation of genetic algorithm
    /// </summary>
    class Optimizer
    {
        public Optimizer()
        {
            rnd = new Random();
        }

        #region GLOBAL VARIABLES
        /// <summary>
        /// Optimization mode (1 - maximalization, 0 - minimalization)
        /// </summary>
        public static int optimizationMode { get; set; }

        /// <summary>
        /// Population size
        /// </summary>
        public static int PopulationSize { get; set; }

        /// <summary>
        /// Number of reproductive cycles
        /// </summary>
        public static int ReproductiveCyclesNumber { get; set; }

        /// <summary>
        /// Probability of mutation
        /// </summary>
        public static float MutationProbability { get; set; }

        /// <summary>
        /// Probability of crossover
        /// </summary>
        public static float CrossoverProbability { get; set; }

        /// <summary>
        /// Size of tournament
        /// </summary>
        public static int TournamentSize { get; set; }

        /// <summary>
        /// Additional stop criterion 
        /// (a number of reproduction cycles ehich not generate improvement)
        /// </summary>
        public static int StopCriterion { get; set; }

        /// <summary>
        /// Population of individuals
        /// </summary>
        private List<List<string>> Population;

        /// <summary>
        /// New population of individuals
        /// </summary>
        private List<List<string>> NewPopulation;

        /// <summary>
        /// Individuals scores
        /// </summary>
        private List<double> PopulationScores;

        /// <summary>
        /// Scores for new population after selection
        /// </summary>
        private List<double> NewPopulationScores;

        /// <summary>
        /// Best score
        /// </summary>
        private double BestScore;

        /// <summary>
        /// Best individual
        /// </summary>
        private List<string> BestIndividual;

        /// <summary>
        /// Codons grouped by amino acids
        /// </summary>
        private Dictionary<string, List<string>> codonGroups;

        /// <summary>
        /// Random declaration
        /// </summary>
        Random rnd;

        System.IO.StreamWriter outSeq;

        #endregion

        #region METHODS
        /// <summary>
        /// Method for grouping codons by amino acid
        /// (key: amino acid, values: list of sysnonymous codons)
        /// </summary>
        private void aminoToCodon()
        {
            codonGroups = new Dictionary<string, List<string>>();
            codonGroups = SeqParser.codonToAmino.GroupBy(x => x.Value)
                .ToDictionary(x => x.Key, x => x.Select(i => i.Key).ToList());
        }

        /// <summary>
        /// Method for codon randomization
        /// Randomizaton of codon for given amino acid
        /// </summary>
        /// <param name="amino"></param>
        /// <returns></returns>
        private string randomizeCodon(string amino)
        {
            string codon = codonGroups[amino][rnd.Next(0, codonGroups[amino].Count())];
            return codon;
        }

        /// <summary>
        /// method for updating best individual and best score
        /// </summary>
        private void updateBestIndividual()
        {
            // for function maximalization
            if (optimizationMode == 1)
            {
                for (int i = 0; i < PopulationScores.Count(); i++)
                {
                    if (PopulationScores[i] > BestScore)
                    {
                        BestScore = PopulationScores[i];
                        BestIndividual.Clear();
                        foreach (string c in Population[i])
                        {
                            BestIndividual.Add(c);
                        }
                    }
                }
            }
            // for function minimalization
            if (optimizationMode == 0)
            {
                for (int i = 0; i < PopulationScores.Count(); i++)
                {
                    if (PopulationScores[i] < BestScore)
                    {
                        BestScore = PopulationScores[i];
                        BestIndividual.Clear();
                        foreach (string c in Population[i])
                        {
                            BestIndividual.Add(c);
                        }
                    }
                }
            }
        }

        /// <summary>
        /// Method for generation of initial population
        /// </summary>
        /// <param name="AminoORFseq"></param>
        private void generateInitialPopulation(List<string> AminoORFseq, string stopCodon, System.IO.StreamWriter outSeq)
        {
            //outSeq.WriteLine("====INITIAL POPULATION====");
            // new population of individuals and scores initialization, new best individual initialization
            Population = new List<List<string>>();
            PopulationScores = new List<double>();
            BestIndividual = new List<string>();

            // temporary variables
            List<string> tempIndividual;
            string tempCodon;

            for (int i = 0; i < PopulationSize; i++)
            {
                // new individual initialization
                tempIndividual = new List<string>();

                // randomization of codons for given amino acid sequence
                foreach (string amino in AminoORFseq)
                {
                    if (amino != "/")
                    {
                        tempCodon = randomizeCodon(amino);
                        tempIndividual.Add(tempCodon);
                    }
                    else
                    {
                        tempIndividual.Add(stopCodon);
                    }
                }

                // adding new individual to population
                Population.Add(tempIndividual);
                // calculation of score for nex individual 
                PopulationScores.Add(ORF.CPBcalculator(tempIndividual));

                /*foreach (string c in tempIndividual)
                {
                    outSeq.Write(c + " ");
                }
                outSeq.WriteLine(PopulationScores.Last());
                outSeq.WriteLine();*/
            }

            BestScore = PopulationScores[0];
            foreach (string c in Population[0])
            {
                BestIndividual.Add(c);
            }
            updateBestIndividual();
        }

        /// <summary>
        /// Method for selection of the best individuals for crossover (tourtnament selection)
        /// </summary>
        private void selectIndividualsForCrossover(int TournamentSize, System.IO.StreamWriter outSeq)
        {
            // new population and new scores initialization
            NewPopulation = new List<List<string>>();
            NewPopulationScores = new List<double>();

            // temporary variables
            int IndividualIdx = 0;
            int BestIndividualIdx = 0;

            for (int i = 0; i < PopulationSize; i++)
            {
                // randomization of individuals for tournament
                for (int j = 0; j < TournamentSize; j++)
                {
                    // randomization of current individual
                    IndividualIdx = rnd.Next(0, PopulationSize);

                    // setting first best score
                    if (j == 0)
                    {
                        BestIndividualIdx = IndividualIdx;
                    }
                    // checking if current individual score is better than current the best
                    else
                    {
                        // optimization mode
                        if (optimizationMode == 1)
                        {
                            if (PopulationScores[IndividualIdx] > PopulationScores[BestIndividualIdx])
                            {
                                BestIndividualIdx = IndividualIdx;
                            }
                        }
                        // deoptimization mode
                        if (optimizationMode == 0)
                        {
                            if (PopulationScores[IndividualIdx] < PopulationScores[BestIndividualIdx])
                            {
                                BestIndividualIdx = IndividualIdx;
                            }
                        }
                    }

                }
                //Console.WriteLine(BestIndividualIdx + " ");

                // adding best individual and score to new population
                NewPopulation.Add(Population[BestIndividualIdx]);
                NewPopulationScores.Add(PopulationScores[BestIndividualIdx]);
            }

            // clearing current population
            Population.Clear();
            PopulationScores.Clear();
            //outSeq.WriteLine();
            //outSeq.WriteLine("====SELECTED POPULATION====");
            // setting new selected population as current population
            for (int i = 0; i < PopulationSize; i++)
            {
                Population.Add(NewPopulation[i]);
                PopulationScores.Add(NewPopulationScores[i]);
                /*foreach (string c in Population.Last())
                {
                    outSeq.Write(c + " ");
                }
                outSeq.WriteLine(PopulationScores.Last());
                outSeq.WriteLine();*/
            }
        }

        /// <summary>
        /// Method for crossover (uniform crossover)
        /// </summary>
        private void crossover(System.IO.StreamWriter outSeq)
        {
            //outSeq.WriteLine();
            //outSeq.WriteLine("====CROSSOVER====");
            // temporary variables
            // first parent and second parent indexes
            int FirstParentIdx, SecondParentIdx;
            // new individuals 
            List<string> FirstNewIndividual, SecondNewIndividual;

            // crossover mask 
            List<int> CrossoverMask;
            int CrossoverMaskSize = Population[0].Count();

            //clearing new population and scores
            NewPopulation.Clear();
            NewPopulationScores.Clear();
            double end;

            if (Math.Round(PopulationSize * CrossoverProbability) % 2 != 0)
            {
                end = PopulationSize - Math.Round(PopulationSize * CrossoverProbability) + 1;
            }
            else
            {
                end = PopulationSize - Math.Round(PopulationSize * CrossoverProbability);
            }

            for (int i = 0; i < end; i++)
            {
                FirstParentIdx = rnd.Next(0, Population.Count());
                NewPopulation.Add(Population[FirstParentIdx]);
                NewPopulationScores.Add(PopulationScores[FirstParentIdx]);
                Population.RemoveAt(FirstParentIdx);
                PopulationScores.RemoveAt(FirstParentIdx);
            }

            // randomization of parents for cross over
            for (int i = 0; i < (PopulationSize - end) / 2; i++)
            {
                outSeq.WriteLine();
                FirstParentIdx = rnd.Next(0, Population.Count());
                SecondParentIdx = rnd.Next(0, Population.Count());
                //outSeq.WriteLine("PARENTS");
                //outSeq.Write("First parent: " + FirstParentIdx + " \n");
                //outSeq.Write("Second parent: " + SecondParentIdx + "\n");

                /*foreach (string c in Population[FirstParentIdx])
                {
                    outSeq.Write(c + " ");
                }

                outSeq.WriteLine();

                foreach (string c in Population[SecondParentIdx])
                {
                    outSeq.Write(c + " ");
                }*/

                //outSeq.WriteLine();

                // rerandomization if parent index was repeated
                while (FirstParentIdx == SecondParentIdx)
                {
                    SecondParentIdx = rnd.Next(0, Population.Count());
                    //outSeq.Write("Second parent again: " + SecondParentIdx + "\n");
                }
                //outSeq.WriteLine();
                // new crossover mask initialization
                CrossoverMask = new List<int>();
                //outSeq.WriteLine("CROSSOVER MASK");
                for (int x = 0; x < CrossoverMaskSize; x++)
                {
                    CrossoverMask.Add(rnd.Next(0, 2));
                    //outSeq.Write(CrossoverMask.Last() + "   ");
                }
                //outSeq.WriteLine();

                // new individuals initialization
                FirstNewIndividual = new List<string>();
                SecondNewIndividual = new List<string>();

                // creation of new individuals using the crossover mask
                for (int x = 0; x < CrossoverMaskSize; x++)
                {
                    if (CrossoverMask[x] == 0)
                    {
                        FirstNewIndividual.Add(Population[FirstParentIdx][x]);
                        SecondNewIndividual.Add(Population[SecondParentIdx][x]);
                    }
                    if (CrossoverMask[x] == 1)
                    {
                        FirstNewIndividual.Add(Population[SecondParentIdx][x]);
                        SecondNewIndividual.Add(Population[FirstParentIdx][x]);
                    }
                }
                //outSeq.WriteLine();
                //outSeq.WriteLine("First new individual");
                /*foreach (string c in FirstNewIndividual)
                {
                    outSeq.Write(c + " ");
                }

                outSeq.WriteLine();
                outSeq.WriteLine("Second new individual");
                foreach (string c in SecondNewIndividual)
                {
                    outSeq.Write(c + " ");
                }*/

                //outSeq.WriteLine();

                // creating new population with new individuals and new scores
                NewPopulation.Add(FirstNewIndividual);
                NewPopulation.Add(SecondNewIndividual);
                NewPopulationScores.Add(ORF.CPBcalculator(FirstNewIndividual));
                NewPopulationScores.Add(ORF.CPBcalculator(SecondNewIndividual));

                // removing "used" parents
                if (FirstParentIdx > SecondParentIdx)
                {
                    Population.RemoveAt(FirstParentIdx);
                    Population.RemoveAt(SecondParentIdx);
                }
                else
                {
                    Population.RemoveAt(SecondParentIdx);
                    Population.RemoveAt(FirstParentIdx);
                }
                //Console.WriteLine("Population: " + Population.Count());
            }

            PopulationScores.Clear();

            for (int i = 0; i < NewPopulation.Count(); i++)
            {
                Population.Add(NewPopulation[i]);
                PopulationScores.Add(NewPopulationScores[i]);
            }

            // updating best individual
            updateBestIndividual();

            //outSeq.WriteLine();
            //outSeq.WriteLine("====AFTER CROSSOVER====");

            /*for (int i = 0; i < PopulationSize; i++)
            {
                foreach (string c in Population[i])
                {
                    outSeq.Write(c + " ");
                }

                outSeq.WriteLine(PopulationScores[i]);
                outSeq.WriteLine();
            }*/
        }

        /// <summary>
        /// Method for mutation
        /// </summary>
        private void mutate(List<string> Individual, int IndividualIdx)
        {
            int codonIdx = rnd.Next(0, Individual.Count()-1);
            string amino = SeqParser.codonToAmino[Individual[codonIdx]];

            // replacing randomized codon and recalculating of score
            Individual[codonIdx] = randomizeCodon(amino);
            PopulationScores[IndividualIdx] = ORF.CPBcalculator(Population[IndividualIdx]);
        }

        /// <summary>
        /// ORF optimization
        /// </summary>
        /// <param name="ORFSeq"></param>
        /// <param name="AminoORFseq"></param>
        /// <param name="optimizationMode"></param>
        /// <returns></returns>
        public List<string> optimizeORF(List<string> ORFSeq, List<string> AminoORFseq, object o, DoWorkEventArgs e)
        {
            string stopCodon = ORFSeq.Last();
            using (outSeq = new System.IO.StreamWriter("D:/ga.txt"))
            {
                // codons grouping to dictionary 
                aminoToCodon();

                // initial population generation
                generateInitialPopulation(AminoORFseq, stopCodon, outSeq);
                
                // reproductive cycles
                for (int i = 0; i < ReproductiveCyclesNumber; i++)
                {
                    Thread.Sleep(1);
                    (o as BackgroundWorker).ReportProgress(100 * i / (ReproductiveCyclesNumber - 1));
                    // selection
                    selectIndividualsForCrossover(TournamentSize, outSeq);

                    // crossover
                    crossover(outSeq);

                    // mutation
                    for (int j = 0; j < Math.Round(Population.Count * MutationProbability); j++)
                    {
                        // individual randomization 
                        int individual = rnd.Next(0, Population.Count());
                        // codon position randomization
                        int codon = rnd.Next(0, Population[individual].Count());
                        // translation to amino acid
                        string amino = SeqParser.codonToAmino[Population[individual][codon]];
                        // mutation of given codon
                        mutate(Population[individual], individual);
                    }
                }

                //outSeq.WriteLine("====LAST POPULATION====");
                // best score for last population

                /*for (int i = 0; i < PopulationSize; i++)
                {
                    foreach (string c in Population[i])
                    {
                        outSeq.Write(c + " ");
                    }
                    outSeq.WriteLine();
                    outSeq.WriteLine(PopulationScores[i]);
                }*/

                // updating best individual
                updateBestIndividual();

                return BestIndividual;
            }
        }
        #endregion
    }
        
}
