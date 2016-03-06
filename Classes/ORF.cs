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

namespace CodonOptimizer.Classes
{
    public class ORF
    {
        public ORF()
        {
            ORFSeq = new List<string>();
            AminoORFseq = new List<string>();
        }

        /// <summary>
        /// ORF list of codons 
        /// </summary>
        public List<string> ORFSeq;

        /// <summary>
        /// ORF list of mino acids
        /// </summary>
        public List<string> AminoORFseq;

        /// <summary>
        /// CPB score
        /// </summary>
        public static double CPB;

        /// <summary>
        /// CPBcalculator method
        /// method for CPB calculating
        /// </summary>
        /// <returns></returns>
        public static double CPBcalculator(List<string> orf)
        {
            CPB = 0;
            int n = 0; 
            foreach (string codon in orf)
            {
                if (n != 0 && codon != "TGA" && codon != "TAA" && codon != "TAG")
                {
                    CPB += CCranking.CCranker.CPS[orf[n-1]+codon];
                }
                n++;
            }

            CPB = CPB / (orf.Count - 1);

            return Math.Round(CPB, 4);
        }
    }
}
