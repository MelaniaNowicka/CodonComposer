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
        #region GLOBAL VARIABLES
        /// <summary>
        /// ORF parser
        /// </summary>
        private static ISequenceParser parser;

        /// <summary>
        /// Codons to Amino acids dictionary
        /// </summary>
        public static Dictionary<string, string> codonToAmino = new Dictionary<string, string>() 
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

        public static Dictionary<string, List<string>> enzymesToSequences = new Dictionary<string, List<string>>() 
        {
            {"AatII", new List<string>{"GACGTC"}},
            {"Acc65I",new List<string>{"GGTACC"}},
            {"AccI",new List<string>{"GTAGAC", "GTCTAC", "GTATAC", "GTCGAC"}},
            {"AcII",new List<string>{"AACGTT"}},
            {"AfeI",new List<string>{"AGCGCT"}},
            {"AflII",new List<string>{"CTTAAG"}},
            {"AgeI",new List<string>{"ACCGGT"}},
            {"ApaI",new List<string>{"GGGCCC"}},
            {"ApaLI",new List<string>{"GTGCAC"}},
            {"ApoI",new List<string>{"AAATTC", "GAATTT", "AAATTT", "GAATTC"}},
            {"AscI",new List<string>{"GGCGCGCC"}},
            {"AseI",new List<string>{"ATTAAT"}},
            {"AsiSI",new List<string>{"GCGATCGC"}},
            {"AvrII",new List<string>{"CCTAGG"}},
            {"BamHI",new List<string>{"GGATCC"}},
            {"BclI",new List<string>{"TGATCA"}},
            {"BglII",new List<string>{"AGATCT"}},
            {"Bme1580I",new List<string>{"GGGCAC", "GTGCCC", "GTGCCC", "GGGCAC"}},
            {"BmtI",new List<string>{"GCTAGC"}},
            {"BsaHI",new List<string>{"GACGCC", "GGCGTC", "GACGTC", "GGCGCC"}},
            {"BsiEI",new List<string>{"CGACCG", "CGGTCG", "CGATCG", "CGGCCG"}},
            {"BsiWI",new List<string>{"CGTACG"}},
            {"BspEI",new List<string>{"TCCGGA"}},
            {"BspHI",new List<string>{"TCATGA"}},
            {"BsrGI",new List<string>{"TGTACA"}},
            {"BssHII",new List<string>{"GCGCGC"}},
            {"BstBI",new List<string>{"TTCGAA"}},
            {"BstZ17I",new List<string>{"GTATAC"}},
            {"BtgI",new List<string>{"CCACGG", "CCGTGG", "CCATGG", "CCGCGG"}},
            {"ClaI",new List<string>{"ATCGAT"}},
            {"DraI",new List<string>{"TTTAAA"}},
            {"EaeI",new List<string>{"CGGCCA", "TGGCCG", "CGGCCG", "TGGCCA"}},   
            {"EagI",new List<string>{"CGGCCG"}},
            {"EcoRI",new List<string>{"GAATTC"}},
            {"EcoRV",new List<string>{"GATATC"}},
            {"FseI",new List<string>{"GGCCGGCC"}},
            {"FspI",new List<string>{"TGCGCA"}},
            {"HaeII",new List<string>{"TGCGCA", "GGCGCT", "AGCGCT", "GGCGCC"}},
            {"HincII",new List<string>{"AGCGCC", "GGCGCT", "AGCGCT", "GGCGCC"}},
            {"HindIII",new List<string>{"AAGCTT"}},
            {"HpaI",new List<string>{"GTTAAC"}},
            {"KasI",new List<string>{"GGCGCC"}},
            {"KpnI",new List<string>{"GGTACC"}},
            {"MfeI",new List<string>{"CAATTG"}},
            {"MluI",new List<string>{"ACGCGT"}},
            {"MscI",new List<string>{"TGGCCA"}},
            {"MspA1I",new List<string>{"CAGCGG", "CCGCTG", "CAGCTG", "CCGCGG"}},
            {"MfeI",new List<string>{"CAATTG"}},
            {"MluI",new List<string>{"ACGCGT"}},
            {"MscI",new List<string>{"TGGCCA"}},
            {"MspA1I",new List<string>{"CAGCGG", "CCGCTG", "CAGCTG", "CCGCGG"}},
            {"NaeI",new List<string>{"GCCGGC"}},
            {"NarI",new List<string>{"GGCGCC"}},   
            {"NcoI",new List<string>{"CCATGG"}},
            {"NdeI",new List<string>{"CATATG"}},
            {"NgoMIV",new List<string>{"GCCGGC"}},
            {"NheI",new List<string>{"GCTAGC"}},
            {"NotI",new List<string>{"GCGGCCGC"}},
            {"NruI",new List<string>{"TCGCGA"}},
            {"NsiI",new List<string>{"ATGCAT"}},
            {"NspI",new List<string>{"ACATGC", "GCATGT", "ACATGT", "GCATGC"}},
            {"PacI",new List<string>{"TTAATTAA"}},
            {"PciI",new List<string>{"ACATGT"}},
            {"PmeI",new List<string>{"GTTTAAAC"}},
            {"PmlI",new List<string>{"CACGTG"}},
            {"PsiI",new List<string>{"TTATAA"}},
            {"PspOMI",new List<string>{"GGGCCC"}},
            {"PstI",new List<string>{"CTGCAG"}},
            {"PvuI",new List<string>{"CGATCG"}},
            {"PvuII",new List<string>{"CAGCTG"}},
            {"SacI",new List<string>{"GAGCTC"}},
            {"SacII",new List<string>{"CCGCGG"}},
            {"SalI",new List<string>{"GTCGAC"}},
            {"SbfI",new List<string>{"CCTGCAGG"}},
            {"ScaI",new List<string>{"AGTACT"}},
            {"SfcI",new List<string>{"CTACAG", "CTGCAG", "CTATAG", "CTGTAG"}},
            {"SfoI",new List<string>{"GGCGCC"}},
            {"SgrAI",new List<string>{"CACCGGCG", "CGCCGGTG", "CACCGGTG", "CGCCGGTCG"}},   
            {"SmaI",new List<string>{"CCCGGG"}},
            {"SmlI",new List<string>{"CTCAAG", "CTCGAG", "CTTAAG", "CTTGAG"}},
            {"SnaBI",new List<string>{"TACGTA"}},
            {"SpeI",new List<string>{"ACTAGT"}},
            {"SphI",new List<string>{"GCATGC"}},
            {"SspI",new List<string>{"AATATT"}},
            {"StuI",new List<string>{"AGGCCT"}},
            {"SwaI",new List<string>{"ATTTAAAT"}},
            {"XbaI",new List<string>{"TCTAGA"}},
            {"XhoI",new List<string>{"CTCGAG"}},
            {"XmaI",new List<string>{"CCCGGG"}}
        };
        #endregion

        #region METHODS
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
                aminos.Add(codonToAmino[l]);
            }
            return aminos;
        }
        #endregion
    }
}
