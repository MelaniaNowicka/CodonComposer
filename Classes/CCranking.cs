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
using System.Windows;

namespace CodonOptimizer.Classes
{
    class CCranking
    {
        public CCranking()
        {

        }

        /// <summary>
        /// ORFeome parser
        /// </summary>
        private ISequenceParser parser;

        /// <summary>
        /// ORFeome from file to list
        /// </summary>
        internal List<ISequence> ORFeome;
        
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

        Dictionary<string, char> CodonsToAmino = new Dictionary<string, char>() 
        {
            // Phenylalanine
            {"TTT",'F'},
            {"TTC",'F'},
            // Leucine
            {"TTA",'L'},
            {"TTG",'L'},
            {"CTT",'L'},
            {"CTC",'L'},
            {"CTA",'L'},
            {"CTG",'L'},
            // Isoleucine
            {"ATT",'I'},
            {"ATC",'I'},
            {"ATA",'I'},
            // Methionine
            {"ATG",'M'},
            // Valine
            {"GTT",'V'},
            {"GTC",'V'},
            {"GTA",'V'},
            {"GTG",'V'},
            // Serine
            {"TCT",'S'},
            {"TCC",'S'},
            {"TCA",'S'},
            {"TCG",'S'},
            // Proline
            {"CCT",'P'},
            {"CCC",'P'},
            {"CCA",'P'},
            {"CCG",'P'},
            // Threonine
            {"ACT",'T'},
            {"ACC",'T'},
            {"ACA",'T'},
            {"ACG",'T'},
            //Alanine
            {"GCT",'A'},
            {"GCC",'A'},
            {"GCA",'A'},
            {"GCG",'A'},
            // Tyrosine
            {"TAT",'Y'},
            {"TAC",'Y'},
            // Histidine
            {"CAT",'H'},
            {"CAC",'H'},
            // Glutamine
            {"CAA",'Q'},
            {"CAG",'Q'},
            // Aparagine
            {"AAT",'N'},
            {"AAC",'N'},
            // Lysine
            {"AAA",'K'},
            {"AAG",'K'},
            // Aspartic acid
            {"GAT",'D'},
            {"GAC",'D'},
            // Glutamic acid
            {"GAA",'E'},
            {"GAG",'E'},
            // Cysteine
            {"TGT",'C'},
            {"TGC",'C'},
            // Tryptophan
            {"TGG",'W'},
            // Arginine
            {"CGT",'R'},
            {"CGC",'R'},
            {"CGA",'R'},
            {"CGG",'R'},
            // Serine
            {"AGT",'S'},
            {"AGC",'S'},
            // Arginine
            {"AGA",'R'},
            {"AGG",'R'},
            // Glycine
            {"GGT",'G'},
            {"GGC",'G'},
            {"GGA",'G'},
            {"GGG",'G'}
        }; 

        /// <summary>
        /// Show information about ORFeome in ORFeomeInfoRichTextBox;
        /// </summary>
        public int readORFeome(string file)
        {
            parser = SequenceParsers.FindParserByFileName(file);
            ORFeome = new List<ISequence>();

            // parsing sequence
            try
            {
                using (parser)
                {
                    this.ORFeome = parser.Parse().ToList();
                }
                parser.Close();
            }
            catch (System.IO.FileFormatException)
            {
                string message = "Something went wrong. Probably you tried to use an improper file. Try again. \nFor more information about using Codon Context Ranking check the \"How to use\" page.";
                ModernDialog.ShowMessage(message.ToString(), "Warning", MessageBoxButton.OK);
            }

            return ORFeome.Count;
        }

        /// <summary>
        /// Returns full fragment sequence as a string. Based on .NET Bio Programming Guide.
        /// </summary>
        /// <returns>Sequence string.</returns>
        private string getString(ISequence seq)
        {
            char[] symbols = new char[seq.Count];
            for (long index = 0; index < seq.Count; index++)
            {
                symbols[index] = (char)seq[index];
            }
            return new String(symbols);
        }

        /// <summary>
        /// Counting codon pairs, codons
        /// </summary>
        private void seqencesToList()
        {
            // temporary sequence
            string seqTemp;
            string aminoPair;
            string amino = "nana";


            // lists initialization
            CodonPairs = new List<string>();
            Codons = new List<string>();
            AminoAcidsPairs = new List<string>();
            AminoAcids = new List<string>();

            // writing to file
            using (System.IO.StreamWriter outSeq = new System.IO.StreamWriter(@"D:/sequences.txt"))
            {
                foreach (ISequence seq in this.ORFeome)
                {
                    seqTemp = getString(seq);

                    // codon pairs substrings
                    for (int i = 0; i < seqTemp.Length - 5; i += 3)
                    {
                        // stop codons elimination
                        if ((seqTemp.Substring(i, 3) != "TGA" && (seqTemp.Substring(i+3, 3)) != "TGA") &&
                            (seqTemp.Substring(i, 3) != "TAA" && (seqTemp.Substring(i + 3, 3)) != "TAA") &&
                            (seqTemp.Substring(i, 3) != "TAG" && (seqTemp.Substring(i + 3, 3)) != "TAG"))
                        {
                            // adding codons pairs
                            //outSeq.WriteLine(seqTemp.Substring(i, 6) + "\n");
                            this.CodonPairs.Add(seqTemp.Substring(i, 6));

                            // adding amino acids pairs
                            aminoPair = this.CodonsToAmino[seqTemp.Substring(i, 3)].ToString() 
                                        + this.CodonsToAmino[seqTemp.Substring(i + 3, 3)].ToString();

                            this.AminoAcidsPairs.Add(aminoPair);
                            //outSeq.WriteLine(aminoPair + "\n");
                        }

                    }

                    // codon substrings
                    for (int i = 0; i < seqTemp.Length - 2; i += 3)
                    {
                        if (seqTemp.Substring(i, 3) != "TGA" &&
                            seqTemp.Substring(i, 3) != "TAA" &&
                            seqTemp.Substring(i, 3) != "TAG")
                        {
                            // adding codons
                            //outSeq.WriteLine(seqTemp.Substring(i, 3) + "\n");
                            this.Codons.Add(seqTemp.Substring(i, 3));

                            //adding amino acids
                            amino = this.CodonsToAmino[seqTemp.Substring(i, 3)].ToString();
                            this.AminoAcids.Add(amino);
                            //outSeq.WriteLine(amino + "\n");
                            
                        }
                    }

                    outSeq.WriteLine("\n\n");
                }
            }
        }

        private void codonsToAminoAcids()
        {
            foreach (string cp in CodonPairs)
            {

            }
        }
        
        /// <summary>
        /// CPS counter
        /// </summary>
        internal void countCPS()
        {

            // sequencesToList method initialization
            seqencesToList();

            // counting elements
            var CodonPairCounts = CodonPairs.GroupBy(i => i);
            var CodonCounts = Codons.GroupBy(i => i);

            var AminoAcidPairCounts = AminoAcidsPairs.GroupBy(i => i);
            var AminoAcidCounts = AminoAcids.GroupBy(i => i);

            using (System.IO.StreamWriter outFile = new System.IO.StreamWriter(@"D:/codons.txt"))
            {
                foreach (var c in CodonPairCounts)
                {
                    outFile.WriteLine(c.Key + " " + c.Count());
                }
            }

            using (System.IO.StreamWriter outFile = new System.IO.StreamWriter(@"D:/aminos.txt"))
            {
                foreach (var c in AminoAcidPairCounts)
                {
                    outFile.WriteLine(c.Key + " " + c.Count());
                }
            }

        }
    }
}
