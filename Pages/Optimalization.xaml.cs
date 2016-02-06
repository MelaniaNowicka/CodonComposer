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

namespace CodonOptimizer.Pages
{
    /// <summary>
    /// Interaction logic for Optimalization.xaml
    /// </summary>
    public partial class Optimalization : UserControl
    {
        public Optimalization()
        {
            ORF = new ORF();
            InitializeComponent();
        }

        /// <summary>
        /// CCranking object
        /// </summary>
        internal static ORF ORF { get; set; }

        /// <summary>
        /// OpenFileDialog object
        /// </summary>
        private Microsoft.Win32.OpenFileDialog openFileDialog;

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
        /// LoadORFButton Click event handler
        /// </summary>
        /// <param name="sender"></param>
        /// <param name="e"></param>
        private void LoadORFButton_Click(object sender, RoutedEventArgs e)
        {
            // openFileDialog method initialization
            initializeOpenFileDialog();

             // show openFileDialog file dialog
            Nullable<bool> openResult = openFileDialog.ShowDialog();

            if (openResult == true)
            {
                string file = openFileDialog.FileName; // file handler

                // parseORF method initialization
                //ORF.parseORF(file);

                BeforeOptimalizationDataGrid.ItemsSource = ORF.ORFseq;
            }
            else
            {
                // modern dialog initialization
                string message = "Something went wrong. Probably you tried to use an improper file. Try again. \nFor more information about using Optimalizator check the \"How to use\" page.";
                ModernDialog.ShowMessage(message.ToString(), "Warning", MessageBoxButton.OK);
            }
        }

    }
}
