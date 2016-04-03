using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Documents;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Navigation;
using System.Windows.Shapes;
using CodonOptimizer.Classes;
using Bio;
using Bio.IO;
using FirstFloor.ModernUI.Windows.Controls;
using System.Threading;
using System.Windows.Forms;
using System.ComponentModel;
using CodonOptimizer.Pages;

namespace CodonOptimizer.Classes
{
    class CCranker
    {
        public CCranker()
        {
            CDScount = 0;
            ORFeome = new List<string>();
            CPS = new Dictionary<string, double>();
        }

        #region GLOBAL VARIABLES
        /// <summary>
        /// ORFeome from file to list
        /// </summary>
        internal List<string> ORFeome;

        /// <summary>
        /// CDS count
        /// </summary>
        internal int CDScount;

        /// <summary>
        /// CodonPairs list object
        /// </summary>
        private List<string> CodonPairs;

        /// <summary>
        /// Codons list object
        /// </summary>
        private List<string> Codons;

        /// <summary>
        /// AminoAcidsPairs object
        /// </summary>
        private List<string> AminoAcidsPairs;

        /// <summary>
        /// AminoAcids object
        /// </summary>
        private List<string> AminoAcids;

        /// <summary>
        /// CPS dictionary
        /// </summary>
        internal Dictionary<string, double> CPS;

        /// <summary>
        /// CC ranking file name
        /// </summary>
        public string fileName;

        /// <summary>
        /// Path for saving files
        /// </summary>
        public string path;

        /// <summary>
        /// Dictionary CodonPairCounts declaration
        /// </summary>
        Dictionary<string, int> CodonPairCounts;

        /// <summary>
        /// Dictionary CodonCounts declaration
        /// </summary>
        Dictionary<string, int> CodonCounts;

        /// <summary>
        /// Dictionary AminoAcidPairCounts declaration
        /// </summary>
        Dictionary<string, int> AminoAcidPairCounts;

        /// <summary>
        /// Dictionary AminoAcidCounts declaration
        /// </summary>
        Dictionary<string, int> AminoAcidCounts;
        #endregion

        #region METHODS
        /// <summary>
        /// SequencesToList method 
        /// </summary>
        private void sequencesToList()
        {
            // temporary variables
            string aminoPair;
            string amino;
            int n = 0;

            // lists initialization
            CodonPairs = new List<string>();
            Codons = new List<string>();
            AminoAcidsPairs = new List<string>();
            AminoAcids = new List<string>();
            int stop = 0;

                foreach (string codon in this.ORFeome)
                {
                    // stop codons elimination
                    if ((codon != "TGA" && codon != "TGA") &&
                        (codon != "TAA" && codon != "TAA") &&
                        (codon != "TAG" && codon != "TAG"))
                    {

                        if (n != 0)
                        {
                            // adding codons pairs
                            if ((ORFeome[n - 1] != "TGA" && ORFeome[n - 1] != "TGA") &&
                                (ORFeome[n - 1] != "TAA" && ORFeome[n - 1] != "TAA") &&
                                (ORFeome[n - 1] != "TAG" && ORFeome[n - 1] != "TAG"))
                            {
                                this.CodonPairs.Add(ORFeome[n - 1] + codon);

                                // adding amino acids pairs
                                aminoPair = SeqParser.codonToAmino[ORFeome[n - 1]].ToString()
                                            + SeqParser.codonToAmino[codon].ToString();

                                this.AminoAcidsPairs.Add(aminoPair);
                            }
                        }
                        // adding codons
                        //outSeq.WriteLine(seqTemp.Substring(i, 3) + "\n");
                        this.Codons.Add(codon);

                        //adding amino acids
                        amino = SeqParser.codonToAmino[codon].ToString();
                        this.AminoAcids.Add(amino);
                    }
                    else { stop++; }
                    n++;
                }
        }

        /// <summary>
        /// elemCounter method
        /// method for counting elements (codons, codon pairs, aminos, amino pairs) in lists
        /// </summary>
        public void elemCounter()
        {

            // sequencesToList method initialization
            sequencesToList();

            // counting elements to dictionaries
            CodonPairCounts = CodonPairs.GroupBy(i => i)
                .ToDictionary(i => i.Key, i => i.Count());

            CodonCounts = Codons.GroupBy(i => i)
                .ToDictionary(i => i.Key, i => i.Count());

            AminoAcidPairCounts = AminoAcidsPairs.GroupBy(i => i)
                .ToDictionary(i => i.Key, i => i.Count());

            AminoAcidCounts = AminoAcids.GroupBy(i => i)
                .ToDictionary(i => i.Key, i => i.Count());

            // codon pairs counting results to file
            using (System.IO.StreamWriter outFile = new System.IO.StreamWriter(path + @"/cpCounts.txt"))
            {
                foreach (var cp in CodonPairCounts)
                {
                    outFile.WriteLine(cp.Key + " " + cp.Value);
                }
            }

            // amino acid pairs counting results to file
            using (System.IO.StreamWriter outFile = new System.IO.StreamWriter(path + @"/apCounts.txt"))
            {
                foreach (var cp in AminoAcidPairCounts)
                {
                    outFile.WriteLine(cp.Key + " " + cp.Value);
                }
            }

            // codons counting results to file
            using (System.IO.StreamWriter outFile = new System.IO.StreamWriter(path + @"/cCounts.txt"))
            {
                foreach (var cp in CodonCounts)
                {
                    outFile.WriteLine(cp.Key + " " + cp.Value);
                }
            }

            // amino acids counting results to file
            using (System.IO.StreamWriter outFile = new System.IO.StreamWriter(path + @"/aaCounts.txt"))
            {
                foreach (var cp in AminoAcidCounts)
                {
                    outFile.WriteLine(cp.Key + " " + cp.Value);
                }
            }

        }

        /// <summary>
        /// CPScounter method
        /// method for calculating CPS
        /// </summary>
        public void CPScalculator(object o, DoWorkEventArgs args)
        {
            // temporary variables
            string aminoPair;
            string codonPair;
            string c1;
            string c2;
            double CPScore;
            int counter = 0;

            // CPS dictionary declaration
            this.CPS = new Dictionary<string, double>();

            // temporary variables
            // frequencies
            // fab - frequency of codon pair
            // fxy - frequency of amino pair
            // fa - frequency of first codon 
            // fb - frequency of second codon
            // fx - frequency of first amino acid
            // fy - frequency of second amino acid
            int fab = 0, fxy = 0, fa = 0, fb = 0, fx = 0, fy = 0;

            foreach (var cp in CodonPairCounts)
            {
                // thread handling
                Thread.Sleep(1);
                (o as BackgroundWorker).ReportProgress(100 * counter / (CodonPairCounts.Count - 1));

                // fab definition 
                fab = cp.Value;

                codonPair = cp.Key.ToString();
                aminoPair = SeqParser.codonToAmino[codonPair.Substring(0, 3)].ToString() + SeqParser.codonToAmino[codonPair.Substring(3, 3)].ToString();

                // fxy definition
                fxy = AminoAcidPairCounts[aminoPair];

                // codons substrings
                c1 = codonPair.Substring(0, 3).ToString();
                c2 = codonPair.Substring(3, 3).ToString();

                // fa, fb definition

                fa = CodonCounts[c1];
                fb = CodonCounts[c2];

                // fx, fy definition

                fx = AminoAcidCounts[SeqParser.codonToAmino[c1].ToString()];
                fy = AminoAcidCounts[SeqParser.codonToAmino[c2].ToString()];

                // CPS counting
                CPScore = Math.Log((double)fab / ((((double)fa * (double)fb) / ((double)fx * (double)fy)) * (double)fxy));
                CPS.Add(codonPair, CPScore);
                counter++;
            }

            // CPS to file CPScores
            using (System.IO.StreamWriter outFile = new System.IO.StreamWriter(path + @"/CPScores.txt"))
            {
                using (System.IO.StreamWriter outFileCSV = new System.IO.StreamWriter(fileName + @".csv"))
                {
                    int n = 1;
                    foreach (KeyValuePair<string, double> cps in CPS)
                    {
                        outFile.WriteLine(n + ". " + cps.Key + " " + cps.Value);
                        outFileCSV.WriteLine(cps.Key + ";" + cps.Value + ";");
                        n++;
                    }
                }
            }
        }

        #endregion
    }
}
