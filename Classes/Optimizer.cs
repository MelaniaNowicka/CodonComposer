using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
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
            MinimalNc = 1.0;
            MaximalNc = 1.0;
        }

        #region GLOBAL VARIABLES
        /// <summary>
        /// Optimization mode (1 - maximalization, 0 - minimalization)
        /// </summary>
        public static int OptimizationMode { get; set; }

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
        /// lysine homopolymers indexes (start, stop) of long A-homopolymers
        /// </summary>
        private List<KeyValuePair<int,int>> lysineIdx;

        /// <summary>
        /// Flag for A-homopolymers removal
        /// </summary>
        public static bool AHomopolymersRemoval { get; set; }

        /// <summary>
        /// Flag for Nc optimization
        /// </summary>
        public static bool NcRestrictions { get; set; }

        /// <summary>
        /// Flag for Nc maintenance
        /// </summary>
        public static bool MaintainOriginalNc { get; set; }

        /// <summary>
        /// Minimal Nc value
        /// </summary>
        public static double MinimalNc { get; set; }

        /// <summary>
        /// Maximal Nc value
        /// </summary>
        public static double MaximalNc { get; set; } 

        private double minimalNc;

        private double maximalNc;

        /// <summary>
        /// Flag for restriction enzyme sites to removal
        /// </summary>
        public static bool RestrEnzymeSitesToRemoval { get; set; }

        /// <summary>
        /// Enzmes Sites to removal 
        /// </summary>
        public static List<string> EnzymeSitesToRemoval { get; set; }

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

        /// <summary>
        /// 
        /// </summary>
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
            for (int i = 0; i < PopulationScores.Count(); i++)
            {
                // for function maximalization
                if (OptimizationMode == 1)
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
                // for function minimalization
                if (OptimizationMode == 0)
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
        /// 
        /// </summary>
        /// <param name="orf"></param>
        private List<KeyValuePair<int, int>> HomopolymersCheck(List<string> AminoORFseq)
        {
            List<KeyValuePair<int, int>> lysineIdx = new List<KeyValuePair<int, int>>();
            int start = 0, stop = 0;

            for (int i = 1; i < AminoORFseq.Count() - 1; i++)
            {
                if(AminoORFseq[i-1] != "K" && AminoORFseq[i] == "K" && AminoORFseq[i+1] == "K")
                {
                    start = i;
                }

                if (AminoORFseq[i - 1] == "K" && AminoORFseq[i] == "K" && AminoORFseq[i + 1] != "K")
                {
                    stop = i;
                }

                if (start != 0 && stop != 0)
                {
                    lysineIdx.Add(new KeyValuePair<int, int>(start, stop));
                    start = 0; stop = 0;
                }
            }

            if (lysineIdx.Count() == 0)
            {
                AHomopolymersRemoval = false;
            }

            return lysineIdx;
        }

        /// <summary>
        /// Method for homopolymers removal
        /// </summary>
        private void HomopolymersRemove(List<string> orf)
        {
            foreach(KeyValuePair<int, int> li in lysineIdx) 
            {
                for (int i = li.Key; i <= li.Value; i++)
                {
                    orf[i] = "AAG";
                }
            }
        }

        /// <summary>
        /// Method for swapping codons
        /// </summary>
        /// <param name="idx"></param>
        /// <param name="orf"></param>
        /// <param name="range"></param>
        private void changeCodon(int idx, List<string> orf, int range)
        {
            int changeIdx = idx + rnd.Next(0, range);
            var matches = SeqParser.codonToAmino
                            .Where(x => x.Value == SeqParser.codonToAmino[orf[changeIdx]] && !x.Key.Equals(orf[changeIdx]))
                            .Select(x => x.Key);

            if (matches.Count() != 0)
            {
                string elem = matches.ElementAt(rnd.Next(0, matches.Count()));        
                orf[changeIdx] = elem;
            }
        }

        /// <summary>
        /// Method for restriction enzyme sites removal
        /// </summary>
        private bool enzymeSitesRemove(List<string> orf, bool allowed, int i)
        {
            string orfStr = string.Join("", orf);
            int idx;
            allowed = true;

            if (EnzymeSitesToRemoval != null)
            {
                foreach (string enzyme in EnzymeSitesToRemoval)
                {
                    foreach (string site in SeqParser.enzymesToSequences[enzyme])
                    {
                        idx = 0;
                        while (idx >= 0)
                        {
                            //searching for enzyme sites within orf
                            idx = orfStr.IndexOf(site, idx);

                            if (idx != -1)
                            {
                                if (idx % 3 == 0)
                                {
                                    if (site.Length == 6)
                                    {
                                        changeCodon(idx / 3, orf, 2);
                                    }
                                    if (site.Length == 8)
                                    {
                                        changeCodon(idx / 3, orf, 3);
                                    }
                                }
                                if (idx % 3 == 1)
                                {
                                    if (site.Length == 6)
                                    {
                                        changeCodon((idx - 1) / 3, orf, 2);
                                    }
                                    if (site.Length == 8)
                                    {
                                        changeCodon((idx - 1) / 3, orf, 3);
                                    }
                                }
                                if (idx % 3 == 2)
                                {
                                    if (site.Length == 6)
                                    {
                                        changeCodon((idx + 1) / 3, orf, 2);
                                    }
                                    if (site.Length == 8)
                                    {
                                        changeCodon((idx + 1) / 3, orf, 3);
                                    }
                                }
                                idx++;
                                allowed = false;
                            }
                        }
                    }
                }
            }
            return allowed;
        }

        /// <summary>
        /// Method for generation of initial population
        /// </summary>
        /// <param name="AminoORFseq"></param>
        private void generateInitialPopulation(ORF orf, string stopCodon, System.IO.StreamWriter outSeq)
        {
            bool allowed;
            // new population of individuals and scores initialization, new best individual initialization
            Population = new List<List<string>>();
            PopulationScores = new List<double>();
            BestIndividual = new List<string>();

            // temporary variables
            List<string> tempIndividual;
            string tempCodon;

            if (MaintainOriginalNc != true)
            {
                for (int i = 0; i < PopulationSize; i++)
                {
                    allowed = false;
                    // new individual initialization
                    tempIndividual = new List<string>();

                    // randomization of codons for given amino acid sequence
                    foreach (string amino in orf.aminoORFseq)
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

                    //Restriction enzymes sites removal
                    if (RestrEnzymeSitesToRemoval == true)
                    {
                        while (allowed == false)
                        {
                            allowed = enzymeSitesRemove(tempIndividual, allowed, i);
                        }
                    }

                    //Homopolymers removal
                    if (AHomopolymersRemoval == true)
                    {
                        HomopolymersRemove(tempIndividual);

                        if (RestrEnzymeSitesToRemoval == true)
                        {
                            allowed = false;

                            while (allowed == false)
                            {
                                allowed = enzymeSitesRemove(tempIndividual, allowed, i);
                            }
                        }
                    }

                    Population.Add(tempIndividual);
                    PopulationScores.Add(ORF.CPBcalculator(tempIndividual));

                    BestScore = PopulationScores[0];
                    foreach (string c in Population[0])
                    {
                        BestIndividual.Add(c);
                    }
                    updateBestIndividual();

                }
            }

            if (MaintainOriginalNc == true)
            {
                for (int i = 0; i < PopulationSize; i++)
                {
                    Population.Add(new List<string>(orf.orfSeq));
                    if (OptimizationMode == 1)
                    {
                        PopulationScores.Add(ORF.MultiScore(orf.orfSeq, orf.aminoAcidCounts, minimalNc, maximalNc, 1));
                    }
                    if (OptimizationMode == 0)
                    {
                        PopulationScores.Add(ORF.MultiScore(orf.orfSeq, orf.aminoAcidCounts, minimalNc, maximalNc, 0));
                    }
                }

                for (int i = 0; i < Population.Count(); i++)
                {
                    for (int j = 0; j < 10; j++)
                    {
                        // individual randomization 
                        int individual = rnd.Next(0, Population.Count());

                        mutate(Population[individual], individual, orf);
                    }
                    
                }
                
                BestScore = PopulationScores[0];
                foreach (string c in Population[0])
                {
                    BestIndividual.Add(c);
                }
                updateBestIndividual();
            }

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
                        if (OptimizationMode == 1)
                        {
                            if (PopulationScores[IndividualIdx] > PopulationScores[BestIndividualIdx])
                            {
                                BestIndividualIdx = IndividualIdx;
                            }
                        }
                        // deoptimization mode
                        if (OptimizationMode == 0)
                        {
                            if (PopulationScores[IndividualIdx] < PopulationScores[BestIndividualIdx])
                            {
                                BestIndividualIdx = IndividualIdx;
                            }
                        }
                    }

                }
                // adding best individual and score to new population
                NewPopulation.Add(Population[BestIndividualIdx]);
                NewPopulationScores.Add(PopulationScores[BestIndividualIdx]);
            }

            // clearing current population
            Population.Clear();
            PopulationScores.Clear();

            // setting new selected population as current population
            for (int i = 0; i < PopulationSize; i++)
            {
                Population.Add(NewPopulation[i]);
                PopulationScores.Add(NewPopulationScores[i]);
            }
        }

        /// <summary>
        /// Method for crossover (uniform crossover)
        /// </summary>
        private bool crossover(System.IO.StreamWriter outSeq, ORF orf, double minimalNc, double maximalNc)
        {
            // temporary variables
            // first parent and second parent indexes
            int FirstParentIdx, SecondParentIdx;
            // new individuals 
            List<string> FirstNewIndividual, SecondNewIndividual;
            int stopCounter = 0;
            bool ncFound = true;

            // crossover mask 
            List<int> CrossoverMask;
            int CrossoverMaskSize = Population[0].Count();

            //clearing new population and scores
            NewPopulation.Clear();
            NewPopulationScores.Clear();
            int end;
            bool allowed;

            if (Math.Round(PopulationSize * CrossoverProbability) % 2 != 0)
            {
                end = PopulationSize - (int)Math.Round(PopulationSize * CrossoverProbability) + 1;
            }
            else
            {
                end = PopulationSize - (int)Math.Round(PopulationSize * CrossoverProbability);
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
                
                FirstParentIdx = rnd.Next(0, Population.Count());
                SecondParentIdx = rnd.Next(0, Population.Count());

                // rerandomization if parent index was repeated
                while (FirstParentIdx == SecondParentIdx)
                {
                    SecondParentIdx = rnd.Next(0, Population.Count());
                }

                // new crossover mask initialization
                CrossoverMask = new List<int>();

                for (int x = 0; x < CrossoverMaskSize; x++)
                {
                    CrossoverMask.Add(rnd.Next(0, 2));
                }

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

                //Restricion enzymes sites removal

                if (RestrEnzymeSitesToRemoval == true)
                {
                    allowed = false;

                    while (allowed == false)
                    {
                        allowed = enzymeSitesRemove(FirstNewIndividual, allowed, i);
                    }

                    allowed = false;

                    while (allowed == false)
                    {
                        allowed = enzymeSitesRemove(SecondNewIndividual, allowed, i);
                    }
                }

                //Homopolymers removal
                if (AHomopolymersRemoval == true)
                {
                    HomopolymersRemove(FirstNewIndividual);
                    HomopolymersRemove(SecondNewIndividual);

                    allowed = false;

                    //Restricion enzymes sites removal
                    if (RestrEnzymeSitesToRemoval == true)
                    {
                        allowed = false;

                        while (allowed == false)
                        {
                            allowed = enzymeSitesRemove(FirstNewIndividual, allowed, i);
                        }

                        allowed = false;

                        while (allowed == false)
                        {
                            allowed = enzymeSitesRemove(SecondNewIndividual, allowed, i);
                        }
                    }
                }
                
                if (MaintainOriginalNc == true)
                {
                    /*double scoreNcFirst = ORF.NcCalculator(FirstNewIndividual, orf.aminoAcidCounts);
                    double scoreNcSecond = ORF.NcCalculator(SecondNewIndividual, orf.aminoAcidCounts);
                    
                    if (scoreNcFirst >= minimalNc && scoreNcFirst <= maximalNc && scoreNcSecond >= minimalNc && scoreNcSecond <= maximalNc)
                    {*/
                        // creating new population with new individuals and new scores
                        NewPopulation.Add(FirstNewIndividual);
                        NewPopulation.Add(SecondNewIndividual);
                        if (OptimizationMode == 1) 
                        {
                            NewPopulationScores.Add(ORF.MultiScore(FirstNewIndividual, orf.aminoAcidCounts, minimalNc, maximalNc, 1));
                            NewPopulationScores.Add(ORF.MultiScore(SecondNewIndividual, orf.aminoAcidCounts, minimalNc, maximalNc, 1));
                        }
                        if (OptimizationMode == 0)
                        {
                            NewPopulationScores.Add(ORF.MultiScore(FirstNewIndividual, orf.aminoAcidCounts, minimalNc, maximalNc, 0));
                            NewPopulationScores.Add(ORF.MultiScore(SecondNewIndividual, orf.aminoAcidCounts, minimalNc, maximalNc, 0));
                        }

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
                        //stopCounter = 0;
                        //Console.WriteLine(stopCounter);
                    /*}
                    else 
                    { 
                        i--; 
                        stopCounter++;
                        if (stopCounter == StopCriterion%2)
                        {
                            for (int j = 0; j < 10; j++)
                            {

                            }
                        }
                        Console.WriteLine(stopCounter +" " + scoreNcFirst + " " + scoreNcSecond);
                        if (stopCounter == StopCriterion) 
                        {
                            Console.WriteLine("cross " + stopCounter + " " + StopCriterion);
                            i = (PopulationSize - end) / 2 - 1; 
                            ncFound = false;
                            break;
                        } 
                    }*/

                }
                else
                {
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
                }
            }

            PopulationScores.Clear();

            for (int j = 0; j < NewPopulation.Count(); j++)
            {
                Population.Add(NewPopulation[j]);
                PopulationScores.Add(NewPopulationScores[j]);
            }

            // updating best individual
            updateBestIndividual();
            return ncFound;
        }

        /// <summary>
        /// Method for mutation
        /// </summary>
        private void mutate(List<string> Individual, int IndividualIdx, ORF orf)
        {
            int codonIdx = rnd.Next(0, Individual.Count()-1);
            string amino = SeqParser.codonToAmino[Individual[codonIdx]];

            // replacing randomized codon and recalculating of score
            Individual[codonIdx] = randomizeCodon(amino);

            //Homopolymers Removal
            if (AHomopolymersRemoval == true)
            {
                HomopolymersRemove(Individual);
            }

            if (MaintainOriginalNc == true)
            {
                if (OptimizationMode == 1)
                {
                    PopulationScores[IndividualIdx] = ORF.MultiScore(Population[IndividualIdx], orf.aminoAcidCounts, minimalNc, maximalNc, 1);
                }
                if (OptimizationMode == 0)
                {
                    PopulationScores[IndividualIdx] = ORF.MultiScore(Population[IndividualIdx], orf.aminoAcidCounts, minimalNc, maximalNc, 0);
                }
            }
            else
            {
                PopulationScores[IndividualIdx] = ORF.CPBcalculator(Population[IndividualIdx]);
            }
            

        }

        /// <summary>
        /// ORF optimization
        /// </summary>
        /// <param name="ORFSeq"></param>
        /// <param name="AminoORFseq"></param>
        /// <param name="optimizationMode"></param>
        /// <returns></returns>
        public List<string> optimizeORF(ORF orf, object o, DoWorkEventArgs e)
        {
            string stopCodon = orf.orfSeq.Last();
            int stopCounter = 0;
            double lastScore = 0;
            
            minimalNc = MinimalNc;
            maximalNc = MaximalNc;

            if (MaintainOriginalNc == true)
            {
                minimalNc = ORF.NcCalculator(orf.orfSeq, orf.aminoAcidCounts) - minimalNc;
                maximalNc = ORF.NcCalculator(orf.orfSeq, orf.aminoAcidCounts) + maximalNc;
            }

            //Homopolymers counting
            if (AHomopolymersRemoval == true)
            {
                lysineIdx = HomopolymersCheck(orf.aminoORFseq);
            }

            // codons grouping to dictionary 
            aminoToCodon();

            // initial population generation
            generateInitialPopulation(orf, stopCodon, outSeq);

            // reproductive cycles
            for (int i = 0; i < ReproductiveCyclesNumber; i++)
            {
                // mutation
                for (int j = 0; j < Math.Round(Population.Count * MutationProbability); j++)
                {
                    // individual randomization 
                    int individual = rnd.Next(0, Population.Count());
                    // mutation of given codon
                    mutate(Population[individual], individual, orf);
                }

                if (CrossoverProbability != 0)
                {
                    // selection
                    selectIndividualsForCrossover(TournamentSize, outSeq);

                    // crossover
                    if (!crossover(outSeq, orf, minimalNc, maximalNc)) { i = ReproductiveCyclesNumber - 1; }

                    if (lastScore != BestScore)
                    {
                        lastScore = BestScore;
                        stopCounter = 0;
                    }
                    else
                    {
                        stopCounter++;
                    }

                    if (stopCounter == StopCriterion)
                    {
                        i = ReproductiveCyclesNumber - 1;
                    }
                }

                Thread.Sleep(1);
                (o as BackgroundWorker).ReportProgress(100 * i / (ReproductiveCyclesNumber - 1));
            }
            // updating best individual
            updateBestIndividual();

            return BestIndividual;

        }
        #endregion
    }
        
}
