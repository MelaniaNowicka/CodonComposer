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

namespace CodonOptimizer.Pages
{
    /// <summary>
    /// 
    /// </summary>
    public partial class CCranking : System.Windows.Controls.UserControl
    {
        public CCranking()
        {
            CCranker = new CCranker();
            InitializeComponent();
        }

        /// <summary>
        /// OpenFileDialog object
        /// </summary>
        private Microsoft.Win32.OpenFileDialog openFileDialog;

        /// <summary>
        /// CCranking object
        /// </summary>
        internal static CCranker CCranker { get; set; }

        /// <summary>
        /// Background worker
        /// </summary>
        internal BackgroundWorker worker = new BackgroundWorker();

        /// <summary>
        /// OpenFileDialog initialization method
        /// </summary>
        private void initializeOpenFileDialog()
        {
            this.openFileDialog = new Microsoft.Win32.OpenFileDialog();
            this.openFileDialog.FileName = ""; // default file name
            this.openFileDialog.Filter = "Fasta files|*.fa;*.fas;*.fasta"; // filter files by extension
            this.openFileDialog.Multiselect = false; // Only one file
            this.openFileDialog.Title = "Open ORFeome file..."; // title text
        }

        /// <summary>
        /// addORFeomeButton click event handler
        /// </summary>
        /// <param name="sender"></param>
        /// <param name="e"></param>
        private void addORFeomeButton_Click(object sender, RoutedEventArgs e)
        {
            // richTextBox cleaning
            ORFeomeInfoRichTextBox.Document.Blocks.Clear();
            
            // openFileDialog method initialization
            initializeOpenFileDialog();

            // show openFileDialog file dialog
            Nullable<bool> openResult = openFileDialog.ShowDialog();

            if (openResult == true)
            {
                string file = openFileDialog.FileName; // file handler

                var tupleTemp = SeqParser.sequenceParser(file);

                CCranker.ORFeome = tupleTemp.Item1;
                CCranker.CDScount = tupleTemp.Item2;

                // adding information to ORFeomeInfoRichTextBox
                if (CCranker.CDScount != 0)
                {
                    ORFeomeInfoRichTextBox.AppendText("Number of CDSs: " + CCranker.CDScount);
                }
            }
            else
            {
                // modern dialog initialization
                string message = "Something went wrong. Probably you tried to use an improper file. Try again. \nFor more information about using Codon Context Ranking check the \"How to use\" page.";
                ModernDialog.ShowMessage(message.ToString(), "Warning", MessageBoxButton.OK);
            }
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
            FolderBrowserDialog folderBrowserDialog = new FolderBrowserDialog();
            DialogResult result = folderBrowserDialog.ShowDialog();

            if (result == DialogResult.OK)
            {
                CPSRichTextBox.AppendText("Please, wait a moment for results...\n\n");

                // setting a path
                CCranker.Path = folderBrowserDialog.SelectedPath;

                // background worker initialization, CPR counter initialization
                worker.DoWork += new DoWorkEventHandler(CCranker.countCPS);
                worker.ProgressChanged += new ProgressChangedEventHandler(worker_ProgressChanged);
                worker.RunWorkerCompleted += new RunWorkerCompletedEventHandler(worker_RunWorkerCompleted);
                worker.WorkerReportsProgress = true;

                CCProgressBar.Value = 0;
                
                // Elements counter initialization
                CCranker.elemCounter();
                CPSRichTextBox.AppendText("Calculating...");
                worker.RunWorkerAsync();
            }
            else
            {
                // modern dialog initialization
                string message = "Something went wrong. It is necessary to select files directory before codon context ranking computing. Try again. \nFor more information about using Codon Context Ranking check the \"How to use\" page.";
                ModernDialog.ShowMessage(message.ToString(), "Warning", MessageBoxButton.OK);
            }
            //int n = 1;

            /*foreach (KeyValuePair<string, double> cps in CCranking.CPS)
            {
                CPSRichTextBox.AppendText(n + ". " + cps.Key + ": " + cps.Value + "\n");
                n++;
            }*/
        }

        /// <summary>
        /// Progress Changed handler
        /// </summary>
        /// <param name="sender"></param>
        /// <param name="e"></param>
        void worker_ProgressChanged(object sender, ProgressChangedEventArgs e)
        {
            CCProgressBar.Value = e.ProgressPercentage;
        }

        /// <summary>
        /// Run Worker Completed handler
        /// </summary>
        /// <param name="sender"></param>
        /// <param name="e"></param>
        void worker_RunWorkerCompleted(object sender, RunWorkerCompletedEventArgs e)
        {
            CPSRichTextBox.AppendText("Calculating is completed. Check files in previously selected directory.\n");
            CPSRichTextBox.AppendText("cCounts - single codon counts\naCounts - single amino acids counts\ncpCounts - codon pair counts\napCounts - amino acid pairs counts\nCPSores - Codon Pair Scores\n");
            CPSRichTextBox.AppendText("For more information, please, check 'How To Use' page\n");
        }

    }
}
