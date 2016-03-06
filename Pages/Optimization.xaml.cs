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
using FirstFloor.ModernUI.Windows.Controls;
using System.Data;
using System.IO;
using Microsoft.VisualBasic.FileIO;

namespace CodonOptimizer.Pages
{
    /// <summary>
    /// Interaction logic for Optimalization.xaml
    /// </summary>
    public partial class Optimization : UserControl
    {
        public Optimization()
        {
            ORF = new ORF();
            Optimizer = new Optimizer();
            InitializeComponent();
        }

        /// <summary>
        /// ORF object
        /// </summary>
        internal ORF ORF { get; set; }

        /// <summary>
        /// Optimizer object
        /// </summary>
        internal Optimizer Optimizer { get; set; }

        /// <summary>
        /// A string of extensions for initializeOpenFileDialog method
        /// </summary>
        string extensions;

        /// <summary>
        /// OpenFileDialog object
        /// </summary>
        private Microsoft.Win32.OpenFileDialog openFileDialog;

        /// <summary>
        /// OpenFileDialog initialization method
        /// </summary>
        private void initializeOpenFileDialog(string ext)
        {
            this.openFileDialog = new Microsoft.Win32.OpenFileDialog();
            this.openFileDialog.FileName = ""; // default file name
            this.openFileDialog.Filter = ext; // filter files by extension
            this.openFileDialog.Multiselect = false; // Only one file
            this.openFileDialog.Title = "Open ORFeome file..."; // title text
        }

        /// <summary>
        /// Creating datagrid view for sequence data
        /// </summary>
        private void sequenceToDataGrid(List<string> DNASeq, List<string> AminoSeq, DataGrid dataGrid, TextBox textBox)
        {
            int count = 1;

            // DataTable declaration
            DataTable Data = new DataTable();

            foreach (var k in AminoSeq)
            {
                Data.Columns.Add(new DataColumn(count.ToString()));
                count++;
            }

            Data.Rows.Add(DNASeq.ToArray());
            Data.Rows.Add(AminoSeq.ToArray());

            // 
            dataGrid.ItemsSource = Data.DefaultView;

            // CPBcalculator method initialization
            textBox.Text = ORF.CPBcalculator(DNASeq).ToString();
        }

        /// <summary>
        /// LoadORFButton Click event handler
        /// </summary>
        /// <param name="sender"></param>
        /// <param name="e"></param>
        private void LoadORFButton_Click(object sender, RoutedEventArgs e)
        {
            // extensions setting
            extensions = "Fasta files|*.fa;*.fas;*.fasta";
            // openFileDialog method initialization
            initializeOpenFileDialog(extensions);

            // show openFileDialog file dialog
            Nullable<bool> openResult = openFileDialog.ShowDialog();

            if (openResult == true)
            {
                string file = openFileDialog.FileName; // file handler

                // sequenceParser method initialization
                var tupleTemp = SeqParser.sequenceParser(file);
                ORF.ORFSeq = tupleTemp.Item1;
                // codonToAminoParser method initialization
                ORF.AminoORFseq = SeqParser.codonToAminoParser(ORF.ORFSeq);

                // sequence to data grid
                sequenceToDataGrid(ORF.ORFSeq, ORF.AminoORFseq, BeforeOptimizationDataGrid, CPBscoreTextBox);
            }
            else
            {
                // modern dialog initialization
                string message = "Something went wrong. Probably you tried to use an improper file. Try again. \nFor more information about using Optimalizator check the \"How to use\" page.";
                ModernDialog.ShowMessage(message.ToString(), "Warning", MessageBoxButton.OK);
            }
        }

        private void UploadRankingButton_Click(object sender, RoutedEventArgs e)
        {
            // extensions setting
            extensions = "CSV files|*.csv";
            // openFileDialog method initialization
            initializeOpenFileDialog(extensions);

            // show openFileDialog file dialog
            Nullable<bool> openResult = openFileDialog.ShowDialog();

            if (openResult == true)
            {
                string file = openFileDialog.FileName; // file handler
                using (TextFieldParser parser = new TextFieldParser(openFileDialog.InitialDirectory + openFileDialog.FileName))
                {
                    CCranking.CCranker = new CCranker();
                    parser.SetDelimiters(new string[] { ";" });

                    while (!parser.EndOfData)
                    {
                        string[] fields = parser.ReadFields();
                        CCranking.CCranker.CPS.Add(fields[0], Convert.ToDouble(fields[1]));

                    }
                }
            }
        }

        private void OptimizeORFButton_Click(object sender, RoutedEventArgs e)
        {
            List<string> OptimalizedORF = Optimizer.optimizeORF(ORF.ORFSeq, ORF.AminoORFseq);
            sequenceToDataGrid(OptimalizedORF, ORF.AminoORFseq, AfterOptimalizationDataGrid, NewCPBscoreTextBox);
        }

    }
}
