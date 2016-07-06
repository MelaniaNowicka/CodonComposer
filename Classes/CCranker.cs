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
    public class CCranker
    {
        public CCranker()
        {
            orfeome = new List<string>();
            cps = new Dictionary<string, double>();
        }

        #region GLOBAL VARIABLES
        /// <summary>
        /// ORFeome from file to list
        /// </summary>
        public List<string> orfeome;

        /// <summary>
        /// CodonPairs list object
        /// </summary>
        private List<string> codonPairs;

        /// <summary>
        /// Codons list object
        /// </summary>
        private List<string> codons;

        /// <summary>
        /// AminoAcidsPairs object
        /// </summary>
        private List<string> aminoAcidsPairs;

        /// <summary>
        /// AminoAcids object
        /// </summary>
        private List<string> aminoAcids;

        /// <summary>
        /// CPS dictionary
        /// </summary>
        public Dictionary<string, double> cps;

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
        Dictionary<string, int> codonPairCounts;

        /// <summary>
        /// Dictionary CodonCounts declaration
        /// </summary>
        Dictionary<string, int> codonCounts;

        /// <summary>
        /// Dictionary AminoAcidPairCounts declaration
        /// </summary>
        Dictionary<string, int> aminoAcidPairCounts;

        /// <summary>
        /// Dictionary AminoAcidCounts declaration
        /// </summary>
        Dictionary<string, int> aminoAcidCounts;

        public static int cdsCount { get; set; }

        #endregion

        #region METHODS
        /// <summary>
        /// SequencesToList method 
        /// </summary>
        private void SequencesToList()
        {
            // temporary variables
            string aminoPair;
            string amino;
            int n = 0;

            // lists initialization
            codonPairs = new List<string>();
            codons = new List<string>();
            aminoAcidsPairs = new List<string>();
            aminoAcids = new List<string>();

            foreach (string codon in this.orfeome)
            {
                // stop codons elimination
                if ((codon != "TGA" && codon != "TGA") &&
                    (codon != "TAA" && codon != "TAA") &&
                    (codon != "TAG" && codon != "TAG"))
                {
                    if (n != 0)
                    {
                        // adding codons pairs
                        if ((orfeome[n - 1] != "TGA" && orfeome[n - 1] != "TGA") &&
                            (orfeome[n - 1] != "TAA" && orfeome[n - 1] != "TAA") &&
                            (orfeome[n - 1] != "TAG" && orfeome[n - 1] != "TAG"))
                        {
                            this.codonPairs.Add(orfeome[n - 1] + codon);

                            // adding amino acids pairs
                            aminoPair = SeqParser.codonToAmino[orfeome[n - 1]].ToString()
                                        + SeqParser.codonToAmino[codon].ToString();

                            this.aminoAcidsPairs.Add(aminoPair);
                        }
                    }
                    // adding codons
                    //outSeq.WriteLine(seqTemp.Substring(i, 3) + "\n");
                    this.codons.Add(codon);

                    //adding amino acids
                    amino = SeqParser.codonToAmino[codon].ToString();
                    this.aminoAcids.Add(amino);
                }
                n++;
            }
        }

        /// <summary>
        /// elemCounter method
        /// method for counting elements (codons, codon pairs, aminos, amino pairs) in lists
        /// </summary>
        public void ElemCounter()
        {

            // sequencesToList method initialization
            SequencesToList();

            // counting elements to dictionaries
            codonPairCounts = codonPairs.GroupBy(i => i)
                .ToDictionary(i => i.Key, i => i.Count());

            codonCounts = codons.GroupBy(i => i)
                .ToDictionary(i => i.Key, i => i.Count());

            aminoAcidPairCounts = aminoAcidsPairs.GroupBy(i => i)
                .ToDictionary(i => i.Key, i => i.Count());

            aminoAcidCounts = aminoAcids.GroupBy(i => i)
                .ToDictionary(i => i.Key, i => i.Count());

            // codon pairs counting results to file
            using (System.IO.StreamWriter outFile = new System.IO.StreamWriter(path + @"/cpCounts.csv"))
            {
                foreach (var cp in codonPairCounts)
                {
                    outFile.WriteLine(cp.Key + ";" + cp.Value + ";");
                }
            }

            // amino acid pairs counting results to file
            using (System.IO.StreamWriter outFile = new System.IO.StreamWriter(path + @"/apCounts.csv"))
            {
                foreach (var cp in aminoAcidPairCounts)
                {
                    outFile.WriteLine(cp.Key + ";" + cp.Value + ";");
                }
            }

            // codons counting results to file
            using (System.IO.StreamWriter outFile = new System.IO.StreamWriter(path + @"/cCounts.csv"))
            {
                foreach (var cp in codonCounts)
                {
                    outFile.WriteLine(cp.Key + ";" + cp.Value + ";");
                }
            }

            // amino acids counting results to file
            using (System.IO.StreamWriter outFile = new System.IO.StreamWriter(path + @"/aaCounts.csv"))
            {
                foreach (var cp in aminoAcidCounts)
                {
                    outFile.WriteLine(cp.Key + ";" + cp.Value + ";");
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
            // Progressbar conter
            int counter = 0;
            // CPS dictionary declaration
            this.cps = new Dictionary<string, double>();

            // temporary variables
            // frequencies
            // fab - frequency of codon pair
            // fxy - frequency of amino pair
            // fa - frequency of first codon 
            // fb - frequency of second codon
            // fx - frequency of first amino acid
            // fy - frequency of second amino acid
            int fab = 0, fxy = 0, fa = 0, fb = 0, fx = 0, fy = 0;

            foreach (var cp in codonPairCounts)
            {
                // thread handling
                Thread.Sleep(1);
                (o as BackgroundWorker).ReportProgress(100 * counter / (codonPairCounts.Count - 1));

                // fab definition 
                fab = cp.Value;

                codonPair = cp.Key.ToString();
                aminoPair = SeqParser.codonToAmino[codonPair.Substring(0, 3)].ToString() + SeqParser.codonToAmino[codonPair.Substring(3, 3)].ToString();

                // fxy definition
                fxy = aminoAcidPairCounts[aminoPair];

                // codons substrings
                c1 = codonPair.Substring(0, 3).ToString();
                c2 = codonPair.Substring(3, 3).ToString();

                // fa, fb definition
                fa = codonCounts[c1];
                fb = codonCounts[c2];

                // fx, fy definition
                fx = aminoAcidCounts[SeqParser.codonToAmino[c1].ToString()];
                fy = aminoAcidCounts[SeqParser.codonToAmino[c2].ToString()];

                // CPS counting
                CPScore = Math.Log((double)fab / ((((double)fa * (double)fb) / ((double)fx * (double)fy)) * (double)fxy));
                cps.Add(codonPair, CPScore);
                counter++;
            }

            // CPS to file CPScores
            using (System.IO.StreamWriter outFile = new System.IO.StreamWriter(path + @"/CPScores.csv"))
            {
                foreach (KeyValuePair<string, double> c in cps)
                {
                    outFile.WriteLine(c.Key + ";" + c.Value + ";");
                }
            }
        }

        #endregion
    }
}
