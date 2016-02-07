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
    class ORF
    {

        public ORF()
        {
            ORFseq = new List<string>();
        }

        /// <summary>
        /// List of ORF
        /// </summary>
        internal List<string> ORFseq;

        /// <summary>
        /// CPB score
        /// </summary>
        internal double CPB;

        internal double CPBcalculator()
        {
            CPB = 0;
            int n = 0; 
            foreach (string codon in ORFseq)
            {
                if (n != 0 && codon != "TGA" && codon != "TAA" && codon != "TAG")
                {
                    CPB += CCranking.CCranker.CPS[ORFseq[n-1]+codon];
                }
                n++;
            }

            CPB = CPB / (ORFseq.Count - 1);

            return Math.Round(CPB, 4);
        }

  
    }
}
