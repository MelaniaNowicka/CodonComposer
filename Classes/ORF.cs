using System;
using System.Collections.Generic;
using System.Windows.Controls;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Bio;
using Bio.IO;
using FirstFloor.ModernUI.Windows.Controls;
using System.Windows;
using CodonOptimizer.Classes;
using CodonOptimizer.Pages;
using System.IO;
using Microsoft.VisualBasic.FileIO;

namespace CodonOptimizer.Classes
{
    public class ORF
    {
        public ORF()
        {
            orfSeq = new List<string>();
            aminoORFseq = new List<string>();
        }

        public ORF(ORF orfToCopy)
        {
            orfSeq = new List<string>(orfToCopy.orfSeq);
            aminoORFseq = new List<string>(orfToCopy.aminoORFseq);
        }

        #region GLOBAL VARIABLES
        /// <summary>
        /// ORF list of codons 
        /// </summary>
        public List<string> orfSeq;

        /// <summary>
        /// ORF list of mino acids
        /// </summary>
        public List<string> aminoORFseq;

        /// <summary>
        /// Counts of amino acids for ORF
        /// </summary>
        public Dictionary<string, int> aminoAcidCounts;

        /// <summary>
        /// CPB score
        /// </summary>
        public static double cpb;

        /// <summary>
        /// NC score
        /// </summary>
        public static double nc;

        public static double ns ;
        public static double nk2, nk3, nk4;

        private static int aaCountAll;
        private static int aaCountInd;
        private static double pseudocountsAll;

        public static double multiScore;

        #endregion

        #region METHODS

        /// <summary>
        /// CPBcalculator method
        /// method for CPB calculation
        /// </summary>
        /// <returns></returns>
        public static double CPBcalculator(List<string> orf)
        {
            cpb = 0;
            int i = 0;

            foreach (string codon in orf)
            {
                if (i != (orf.Count - 1) && codon == "TGA" || codon == "TAA" || codon == "TAG")
                {
                    // add an exception
                }

                if (i != 0 && codon != "TGA" && codon != "TAA" && codon != "TAG")
                {
                    bool elemExists = CCranking.CCranker.cps.ContainsKey(orf[i - 1] + codon);
                    if (elemExists == true)
                    {
                        cpb += CCranking.CCranker.cps[orf[i - 1] + codon];
                    }
                    else
                    {
                        cpb += 0;
                    }
                }
                i++;
            }

            cpb = cpb / (double)(orf.Count - 1);

            return Math.Round(cpb, 4);
        }

        public static double MultiScore(List<string> orf, Dictionary<string, int> aminoAcidCounts, double minimalNc, double maximalNc, int optimizationMode)
        {
            double ncScore = NcCalculator(orf, aminoAcidCounts);

            if (optimizationMode == 1)
            {
                if (ncScore < minimalNc)
                {
                    multiScore = CPBcalculator(orf) - (minimalNc - ncScore) * 0.1;
                }

                if (ncScore > maximalNc)
                {
                    multiScore = CPBcalculator(orf) - (ncScore - maximalNc) * 0.1;
                }

                if (ncScore >= minimalNc && ncScore <= maximalNc)
                {
                    multiScore = CPBcalculator(orf);
                }
            }
            if (optimizationMode == 0)
            {
                if (ncScore < minimalNc)
                {
                    multiScore = CPBcalculator(orf) + (minimalNc - ncScore) * 0.1;
                }

                if (ncScore > maximalNc)
                {
                    multiScore = CPBcalculator(orf) + (ncScore - maximalNc) * 0.1;
                }

                if (ncScore >= minimalNc && ncScore <= maximalNc)
                {
                    multiScore = CPBcalculator(orf);
                }
            }

            return Math.Round(multiScore, 4);
        }

        /// <summary>
        /// NCcalculator method
        /// Method for Nc calculation
        /// </summary>
        /// <param name="orf"></param>
        /// <returns></returns>
        public static double NcCalculator(List<string> orf, Dictionary<string, int> aminoAcidCounts)
        {
            nc = 0;
            ns = 0;
            nk2 = 0;
            nk3 = 0;
            nk4 = 0;

            // codons counting
            var codonCounts = orf.GroupBy(i => i).ToDictionary(i => i.Key, i => i.Count());

            ns = GeneticCode.oneFoldFamilies.Count();
            nk2 = NcKiCodonFamiliesCalculator(2, GeneticCode.twoFoldFamilies, aminoAcidCounts, codonCounts);
            nk3 = NcKiCodonFamiliesCalculator(3, GeneticCode.threeFoldFamilies, aminoAcidCounts, codonCounts);
            nk4 = NcKiCodonFamiliesCalculator(4, GeneticCode.fourFoldFamilies, aminoAcidCounts, codonCounts);

            nc = ns + nk2 + nk3 + nk4;
            
            return Math.Round(nc,4);
        }

        /// <summary>
        /// Method for partial NC calculations
        /// </summary>
        /// <returns></returns>
        private static double NcKiCodonFamiliesCalculator(int Ki, Dictionary<string, List<string>> xFoldFamilies, Dictionary<string, int> aminoAcidCounts, Dictionary<string, int> codonCounts)
        {
            aaCountAll = 0;
            aaCountInd = 0;
            pseudocountsAll = 0;

            foreach (KeyValuePair<string, List<string>> family in xFoldFamilies)
            {
                aaCountInd = 0;

                if (aminoAcidCounts.ContainsKey(family.Key))
                {
                    foreach (string codon in family.Value)
                    {
                        if (codonCounts.ContainsKey(codon))
                        {
                            aaCountAll += codonCounts[codon];
                            aaCountInd += codonCounts[codon];
                        }
                    }

                    pseudocountsAll += aaCountInd * FcfCalculator(family.Value, Ki, codonCounts, aaCountInd);
                }
                else
                {
                    pseudocountsAll += 1.0 * FcfCalculator(family.Value, Ki, codonCounts, 1);
                }
            }
            return xFoldFamilies.Count() * aaCountAll / pseudocountsAll;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <returns></returns>
        private static double FcfCalculator(List<string> codonFamily, int codonFamilySize, Dictionary<string, int> codonCounts, int aminoAcidCount)
        {
            double fcf = 0;

            for (int i = 0; i < codonFamily.Count();  i++)
            {
                if (codonCounts.ContainsKey(codonFamily[i]))
                {
                    fcf += Math.Pow((double)(codonCounts[codonFamily[i]] + 1.0) / (double)(aminoAcidCount + codonFamilySize), 2);
                }
                else
                {
                    fcf += Math.Pow(1.0/((double)aminoAcidCount + (double)codonFamilySize), 2);
                }
            }

            return fcf;
        }

        public static bool lysinesHomopolymersCheck(List<string> aminoSeq)
        {
            bool lysinesHomopolymersFound = false;
            int idx;

            idx = 0;
            while (idx >= 0)
            {
                idx = aminoSeq.IndexOf("L", idx);
            }

            return lysinesHomopolymersFound;
        }
        #endregion

    }
}
