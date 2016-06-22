using FirstFloor.ModernUI.Windows.Controls;
using Microsoft.VisualBasic.FileIO;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;

namespace CodonOptimizer.Classes
{
    class GeneticCode
    {
        /// <summary>
        /// Dictionary for 1-fold codon families
        /// </summary>
        public static Dictionary<string, List<string>> oneFoldFamilies;

        /// <summary>
        /// Dictionary for 2-fold codon families
        /// </summary>
        public static Dictionary<string, List<string>> twoFoldFamilies;

        /// <summary>
        /// Dictionary for 3-fold codon families
        /// </summary>
        public static Dictionary<string, List<string>> threeFoldFamilies;

        /// <summary>
        /// Dictionary for 4-fold codon families
        /// </summary>
        public static Dictionary<string, List<string>> fourFoldFamilies;

        public static List<string> stopCodons;

        /// <summary>
        /// Method for genetic code upload
        /// </summary>
        public static void UploadGeneticCode()
        {
            string path = System.IO.Path.Combine(Directory.GetCurrentDirectory(), "gencode_Standard.csv");

            oneFoldFamilies = new Dictionary<string, List<string>>();
            twoFoldFamilies = new Dictionary<string, List<string>>();
            threeFoldFamilies = new Dictionary<string, List<string>>();
            fourFoldFamilies = new Dictionary<string, List<string>>();
            stopCodons = new List<string>();

            using (TextFieldParser parser = new TextFieldParser(path))
            {
                parser.SetDelimiters(new string[] { ";" });

                if (parser.EndOfData)
                {
                    // modern dialog initialization
                    string message = "The genetic code file is empty.";
                    ModernDialog.ShowMessage(message.ToString(), "Warning", MessageBoxButton.OK);
                }
                else
                {
                    while (!parser.EndOfData)
                    {
                        string[] fields = parser.ReadFields();
                        List<string> tmpList = new List<string>();

                        for (int i = 1; i < fields.Count(); i++)
                        {
                            tmpList.Add(fields[i]);
                        }

                        if (tmpList.Count() == 1)
                        {
                            oneFoldFamilies.Add(fields[0], tmpList);
                        }

                        if (tmpList.Count() == 2)
                        {
                            twoFoldFamilies.Add(fields[0], tmpList);
                        }

                        if (tmpList.Count() == 3)
                        {
                            if (fields[0] != "STOP")
                            {
                                threeFoldFamilies.Add(fields[0], tmpList);
                            }
                            else
                            {
                                stopCodons = tmpList;
                            }
                        }

                        if (tmpList.Count() == 4)
                        {
                            fourFoldFamilies.Add(fields[0], tmpList);
                        }
                    }
                }
            }
        }
    }
}
