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
    /// Interaction logic for CCrank.xaml
    /// </summary>
    public partial class CCrank : UserControl
    {
        public CCrank()
        {
            CCranking = new CCranking();

            InitializeComponent();
        }

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
            this.openFileDialog.FileName = ""; // Default file name
            this.openFileDialog.Filter = "Fasta files|*.fa;*.fas;*.fasta"; // Filter files by extension
            this.openFileDialog.Multiselect = false;
            this.openFileDialog.Title = "Open ORFeome file...";
        }

        /// <summary>
        /// CCranking object
        /// </summary>
        internal static CCranking CCranking { get; set; }

        /// <summary>
        /// addORFeomeButton click event handler
        /// </summary>
        /// <param name="sender"></param>
        /// <param name="e"></param>
        private void addORFeomeButton_Click(object sender, RoutedEventArgs e)
        {
            // OpenFileDialog method initialization
            initializeOpenFileDialog();

            // Show openFileDialog file dialog
            Nullable<bool> openResult = openFileDialog.ShowDialog();

            if(openResult == true)
            {
                string file = openFileDialog.FileName;
                int count = 0;
                count = CCranking.readORFeome(file);
                // adding information to ORFeomeInfoRichTextBox
                if (CCranking.ORFeome.Count != 0)
                {
                    ORFeomeInfoRichTextBox.AppendText("Number of CDSs: " + count);
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
            // CPR counter initialization
            CCranking.countCPS();
        }
        
    }
}
