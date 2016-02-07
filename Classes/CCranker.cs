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
        }

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
        /// Path for saving files
        /// </summary>
        internal string Path;

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

        /// <summary>
        /// Counting codon pairs, codons
        /// </summary>
        private void seqencesToList()
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

            // writing to file
            using (System.IO.StreamWriter outSeq = new System.IO.StreamWriter(Path + @"/sequences.txt"))
            {
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
                                aminoPair = SeqParser.CodonsToAmino[ORFeome[n - 1]].ToString()
                                            + SeqParser.CodonsToAmino[codon].ToString();

                                this.AminoAcidsPairs.Add(aminoPair);
                            }
                        }
                        // adding codons
                        //outSeq.WriteLine(seqTemp.Substring(i, 3) + "\n");
                        this.Codons.Add(codon);

                        //adding amino acids
                        amino = SeqParser.CodonsToAmino[codon].ToString();
                        this.AminoAcids.Add(amino);
                    }
                    n++;
                }
            }
        }

        internal void elemCounter()
        {

            // sequencesToList method initialization
            seqencesToList();

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
            using (System.IO.StreamWriter outFile = new System.IO.StreamWriter(Path + @"/cpCounts.txt"))
            {
                foreach (var cp in CodonPairCounts)
                {
                    outFile.WriteLine(cp.Key + " " + cp.Value);
                }
            }

            // amino acid pairs counting results to file
            using (System.IO.StreamWriter outFile = new System.IO.StreamWriter(Path + @"/apCounts.txt"))
            {
                foreach (var cp in AminoAcidPairCounts)
                {
                    outFile.WriteLine(cp.Key + " " + cp.Value);
                }
            }

            // codons counting results to file
            using (System.IO.StreamWriter outFile = new System.IO.StreamWriter(Path + @"/cCounts.txt"))
            {
                foreach (var cp in CodonCounts)
                {
                    outFile.WriteLine(cp.Key + " " + cp.Value);
                }
            }

            // amino acids counting results to file
            using (System.IO.StreamWriter outFile = new System.IO.StreamWriter(Path + @"/aaCounts.txt"))
            {
                foreach (var cp in AminoAcidCounts)
                {
                    outFile.WriteLine(cp.Key + " " + cp.Value);
                }
            }

        }


        /// <summary>
        /// CPS counter
        /// </summary>
        internal void countCPS(object o, DoWorkEventArgs args)
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

            using (System.IO.StreamWriter outFile = new System.IO.StreamWriter(Path + @"/CPSvariables.txt"))
            {
                foreach (var cp in CodonPairCounts)
                {
                    Thread.Sleep(10);
                    (o as BackgroundWorker).ReportProgress(100 * counter / (CodonPairCounts.Count - 1));

                    // writing fab to file
                    //outFile.WriteLine(cp.Key + " " + cp.Value);

                    // fab definition 
                    fab = cp.Value;

                    codonPair = cp.Key.ToString();
                    aminoPair = SeqParser.CodonsToAmino[codonPair.Substring(0, 3)].ToString() + SeqParser.CodonsToAmino[codonPair.Substring(3, 3)].ToString();

                    // fxy definition
                    fxy = AminoAcidPairCounts[aminoPair];

                    // codons substrings
                    c1 = codonPair.Substring(0, 3).ToString();
                    c2 = codonPair.Substring(3, 3).ToString();

                    // fa, fb definition

                    fa = CodonCounts[c1];
                    fb = CodonCounts[c2];

                    // fx, fy definition

                    fx = AminoAcidCounts[SeqParser.CodonsToAmino[c1].ToString()];
                    fy = AminoAcidCounts[SeqParser.CodonsToAmino[c2].ToString()];

                    // writing fxy, fx, fy, fa, fb fo file 
                    //outFile.WriteLine("fxy" + fxy);
                    //outFile.WriteLine("fx" + fx);
                    //outFile.WriteLine("fa" + fa);
                    //outFile.WriteLine("fy" + fy);
                    //outFile.WriteLine("fb" + fb + "\n");

                    // CPS counting
                    CPScore = Math.Log((double)fab / ((((double)fa * (double)fb) / ((double)fx * (double)fy)) * (double)fxy));
                    CPS.Add(codonPair, CPScore);
                    counter++;
                }
            }

            // CPS to file CPScores
            using (System.IO.StreamWriter outFile = new System.IO.StreamWriter(Path + @"/CPScores.txt"))
            {
                int n = 1;

                foreach (KeyValuePair<string, double> cps in CPS)
                {
                    outFile.WriteLine(n + ". " + cps.Key + " " + cps.Value);
                    n++;
                }
            }
        }
    }
}
