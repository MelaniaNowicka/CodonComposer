using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Documents;
using System.Windows.Forms;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Navigation;
using System.Windows.Shapes;
using System.Threading;
using CodonOptimizer.Classes;
using FirstFloor.ModernUI.Windows.Controls;
using System.IO;

namespace CodonOptimizer.Pages
{
    /// <summary>
    /// 
    /// </summary>
    public partial class CCranking : System.Windows.Controls.UserControl
    {
        public CCranking()
        {
            InitializeComponent();
        }

        #region GLOBAL VARIABLES
        /// <summary>
        /// OpenFileDialog object
        /// </summary>
        private Microsoft.Win32.OpenFileDialog openFileDialog;

        /// <summary>
        /// SaveFileDialog object
        /// </summary>
        private Microsoft.Win32.SaveFileDialog saveFileDialog;

        /// <summary>
        /// CCranking object
        /// </summary>
        public static CCranker CCranker { get; set; }

        /// <summary>
        /// Background worker
        /// </summary>
        public BackgroundWorker worker;

        /// <summary>
        /// Flag for checkbox enabling
        /// </summary>
        public static bool currentRankingCheckBoxIsEnabled = false;

        #endregion

        #region METHODS
        /// <summary>
        /// openFileDialog initialization method
        /// </summary>
        private void OpenFileDialogInitialize()
        {
            this.openFileDialog = new Microsoft.Win32.OpenFileDialog();
            this.openFileDialog.FileName = ""; // default file name
            this.openFileDialog.Filter = "Fasta files|*.fa;*.fas;*.fasta"; // filter files by extension
            this.openFileDialog.Multiselect = false; // Only one file
            this.openFileDialog.Title = "Open ORFeome file..."; // title text
        }

        /// <summary>
        /// saveFileDialog
        /// </summary>
        private void SaveFileDialogInitialize()
        {
            this.saveFileDialog = new Microsoft.Win32.SaveFileDialog();
            this.saveFileDialog.FileName = ""; // default file name
            this.openFileDialog.Filter = ".csv files|*.csv"; // filter files by extension
            this.openFileDialog.Multiselect = false; // Only one file
            this.openFileDialog.Title = "Save ranking file..."; // title text
        }

        /// <summary>
        /// addORFeomeButton click event handler
        /// </summary>
        /// <param name="sender"></param>
        /// <param name="e"></param>
        private void AddORFeomeButton_Click(object sender, RoutedEventArgs e)
        {
            // richTextBox cleaning
            ORFeomeInfoRichTextBox.Document.Blocks.Clear();
            CCranker = new CCranker();

            // openFileDialog method initialization
            OpenFileDialogInitialize();

            // show openFileDialog file dialog
            Nullable<bool> openResult = openFileDialog.ShowDialog();

            //if (openResult == true)
            //{
                string file = openFileDialog.FileName; // file handler

                // sequenceParser method initialization
                var tupleTemp = SeqParser.sequenceParser(file);

                CCranker.orfeome = tupleTemp.Item1;
                CCranker.cdsCount = tupleTemp.Item2;

                // adding information to ORFeomeInfoRichTextBox
                if (CCranker.cdsCount != 0)
                {
                    ORFeomeInfoRichTextBox.AppendText(CCranker.cdsCount.ToString());
                }
            /*}
            else
            {
                // modern dialog initialization
                string message = "Something went wrong. Probably you tried to use an improper file. Try again. \nFor more information about using Codon Context Ranking check the \"How to use\" page.";
                ModernDialog.ShowMessage(message.ToString(), "Warning", MessageBoxButton.OK);
            }*/
        }


        /// <summary>
        /// CCrankButton click event handler
        /// </summary>
        /// <param name="sender"></param>
        /// <param name="e"></param>
        private void CCrankButton_Click(object sender, RoutedEventArgs e)
        {
            // richTextBox cleaning
            CPSRichTextBox.Document.Blocks.Clear();
            // directory setting
            SaveFileDialogInitialize();
            Nullable<bool> result = saveFileDialog.ShowDialog();

            if (result == true)
            {
                CPSRichTextBox.AppendText("Please, wait a moment for the results...\n\n");
                worker = new BackgroundWorker();

                // setting a path
                CCranker.fileName = saveFileDialog.FileName;
                CCranker.path = System.IO.Path.Combine(System.IO.Path.GetDirectoryName(saveFileDialog.FileName), CCranker.fileName);
                Directory.CreateDirectory(System.IO.Path.Combine(CCranker.path, CCranker.fileName));

                // Background worker settings
                // background worker initialization, CPR counter initialization
                worker.DoWork += new DoWorkEventHandler(CCranker.CPScalculator);
                worker.ProgressChanged += new ProgressChangedEventHandler(worker_ProgressChanged);
                worker.RunWorkerCompleted += new RunWorkerCompletedEventHandler(worker_RunWorkerCompleted);
                worker.WorkerReportsProgress = true;

                CCProgressBar.Value = 0;
                
                // Elements counter initialization
                CCranker.ElemCounter();
                CPSRichTextBox.AppendText("Calculating...\n");
                worker.RunWorkerAsync();
                currentRankingCheckBoxIsEnabled = true;
            }
            else
            {
                // modern dialog initialization
                string message = "Something went wrong. It is necessary to select files directory before codon context ranking computing. Try again. \nFor more information about using Codon Context Ranking check the \"How to use\" page.";
                ModernDialog.ShowMessage(message.ToString(), "Warning", MessageBoxButton.OK);
            }
        }

        /// <summary>
        /// ProgressChanged handler
        /// </summary>
        /// <param name="sender"></param>
        /// <param name="e"></param>
        void worker_ProgressChanged(object sender, ProgressChangedEventArgs e)
        {
            CCProgressBar.Value = e.ProgressPercentage;
        }

        /// <summary>
        /// RunWorkerCompleted handler
        /// </summary>
        /// <param name="sender"></param>
        /// <param name="e"></param>
        void worker_RunWorkerCompleted(object sender, RunWorkerCompletedEventArgs e)
        {
            CPSRichTextBox.AppendText("Calculation is completed. Check files in previously selected directory.\n");
            CPSRichTextBox.AppendText("cCounts - single codon counts\naCounts - single amino acids counts\ncpCounts - codon pair counts\napCounts - amino acid pairs counts\nCPScores - Codon Pair Scores\n");
            CPSRichTextBox.AppendText("For more information, please, check 'How To Use' page\n");

            // verification of ranking correctness
            // 3721 codon pairs - all possibilities of codon pairs permutations
            if (CCranker.cps.Count() != 3721)
            {
                string message = "Not all possible codon pairs found within the orfeome. Further analysis may be impossible.";
                ModernDialog.ShowMessage(message.ToString(), "Information", MessageBoxButton.OK);
            }
        }
        #endregion
    }
}
