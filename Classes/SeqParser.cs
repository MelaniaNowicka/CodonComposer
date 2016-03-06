using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Bio;
using Bio.IO;
using FirstFloor.ModernUI.Windows.Controls;
using System.Windows;

namespace CodonOptimizer.Classes
{
    public static class SeqParser
    {
        /// <summary>
        /// ORF parser
        /// </summary>
        private static ISequenceParser parser;

        /// <summary>
        /// Codons to Amino acids dictionary
        /// </summary>
        public static Dictionary<string, string> CodonToAmino = new Dictionary<string, string>() 
        {
            // Phenylalanine
            {"TTT","F"},
            {"TTC","F"},
            // Leucine
            {"TTA","L"},
            {"TTG","L"},
            {"CTT","L"},
            {"CTC","L"},
            {"CTA","L"},
            {"CTG","L"},
            // Isoleucine
            {"ATT","I"},
            {"ATC","I"},
            {"ATA","I"},
            // Methionine
            {"ATG","M"},
            // Valine
            {"GTT","V"},
            {"GTC","V"},
            {"GTA","V"},
            {"GTG","V"},
            // Serine
            {"TCT","S"},
            {"TCC","S"},
            {"TCA","S"},
            {"TCG","S"},
            // Proline
            {"CCT","P"},
            {"CCC","P"},
            {"CCA","P"},
            {"CCG","P"},
            // Threonine
            {"ACT","T"},
            {"ACC","T"},
            {"ACA","T"},
            {"ACG","T"},
            //Alanine
            {"GCT","A"},
            {"GCC","A"},
            {"GCA","A"},
            {"GCG","A"},
            // Tyrosine
            {"TAT","Y"},
            {"TAC","Y"},
            // Histidine
            {"CAT","H"},
            {"CAC","H"},
            // Glutamine
            {"CAA","Q"},
            {"CAG","Q"},
            // Aparagine
            {"AAT","N"},
            {"AAC","N"},
            // Lysine
            {"AAA","K"},
            {"AAG","K"},
            // Aspartic acid
            {"GAT","D"},
            {"GAC","D"},
            // Glutamic acid
            {"GAA","E"},
            {"GAG","E"},
            // Cysteine
            {"TGT","C"},
            {"TGC","C"},
            // Tryptophan
            {"TGG","W"},
            // Arginine
            {"CGT","R"},
            {"CGC","R"},
            {"CGA","R"},
            {"CGG","R"},
            // Serine
            {"AGT","S"},
            {"AGC","S"},
            // Arginine
            {"AGA","R"},
            {"AGG","R"},
            // Glycine
            {"GGT","G"},
            {"GGC","G"},
            {"GGA","G"},
            {"GGG","G"},
            // STOP CODONS
            {"TGA", "/"},
            {"TAA", "/"},
            {"TAG", "/"}
        };

        /// <summary>
        /// getString method
        /// Returns full fragment sequence as a string. Based on .NET Bio Programming Guide.
        /// </summary>
        /// <returns>Sequence string.</returns>
        private static string getString(ISequence seq)
        {
            char[] symbols = new char[seq.Count];
            for (long index = 0; index < seq.Count; index++)
            {
                symbols[index] = (char)seq[index];
            }
            return new String(symbols);
        }

        /// <summary>
        /// sequenceParser method
        /// method for parsing sequences from file
        /// return tuple(list of codons, number of cds's)
        /// </summary>
        /// <param name="file"></param>
        public static Tuple<List<string>, int> sequenceParser(string file)
        {

            parser = SequenceParsers.FindParserByFileName(file);
            List<ISequence> sequences = new List<ISequence>();
            List<string> list = new List<string>();

            // temp variables
            string seqTemp;

            // parsing sequence
            try
            {
                using (parser)
                {
                    sequences = parser.Parse().ToList();
                    foreach (ISequence seq in sequences)
                    {
                        // getString method initialization
                        seqTemp = getString(seq);

                        // adding codon substrings
                        for (int i = 0; i < seqTemp.Length - 2; i += 3)
                        {
                            list.Add(seqTemp.Substring(i, 3).ToUpper());
                        }
                    }
                }
                parser.Close();
            }
            catch (System.IO.FileFormatException)
            {
                string message = "Something went wrong. Probably you tried to use an improper file. Try again. \nFor more information about using Codon Context Ranking check the \"How to use\" page.";
                ModernDialog.ShowMessage(message.ToString(), "Warning", MessageBoxButton.OK);
            }
            return new Tuple<List<string>, int>(list, sequences.Count);
        }

        /// <summary>
        /// codonToAminoParser method
        /// method for translating codons to amino acids
        /// </summary>
        /// <param name="list"></param>
        /// <returns></returns>
        public static List<string> codonToAminoParser(List<string> list)
        {
            List<string> aminos = new List<string>();

            // translating codons to amino acids using CodonToAmino dictionary
            foreach (string l in list)
            {
                aminos.Add(CodonToAmino[l]);
            }
            return aminos;
        }
    }
}
