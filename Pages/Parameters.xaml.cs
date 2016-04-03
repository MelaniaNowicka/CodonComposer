using CodonOptimizer.Classes;
using FirstFloor.ModernUI.Windows.Controls;
using System;
using System.Collections.Generic;
using System.Globalization;
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

namespace CodonOptimizer.Pages
{
    /// <summary>
    /// Interaction logic for Parameters.xaml
    /// </summary>
    public partial class Parameters : UserControl
    {
        public Parameters()
        {
            InitializeComponent();
            TournamentSizeTextBox.IsEnabled = false;
            StopCriterionTextBox.IsEnabled = false;
        }

        public Dictionary<string, List<string>> enzymesToSequences = new Dictionary<string, List<string>>() 
        {
            {"AatII", new List<string>{"GACGTC"}},
            {"Acc65I",new List<string>{"GGTACC"}},
            {"AccI",new List<string>{"GTAGAC", "GTCTAC", "GTATAC", "GTCGAC"}},
            {"AclI",new List<string>{"AACGTT"}},
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
            {"MluI",new List<string>{"ACGCGT"}},
            {"MspA1I",new List<string>{"CAGCGG", "CCGCTG", "CAGCTG", "CCGCGG"}},
            {"MfeI",new List<string>{"CAATTG"}},
            {"MscI",new List<string>{"TGGCCA"}},
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

        #region METHODS
        /// <summary>
        /// AdditionalParametersCheckbox_Checked handler
        /// </summary>
        /// <param name="sender"></param>
        /// <param name="e"></param>
        private void AdditionalParametersCheckbox_Checked(object sender, RoutedEventArgs e)
        {
            TournamentSizeTextBox.IsEnabled = true;
            StopCriterionTextBox.IsEnabled = true;
        }

        /// <summary>
        /// AdditionalParametersCheckbox_Unchecked handler
        /// </summary>
        /// <param name="sender"></param>
        /// <param name="e"></param>
        private void AdditionalParametersCheckbox_Unchecked(object sender, RoutedEventArgs e)
        {
            TournamentSizeTextBox.IsEnabled = false;
        }

        /// <summary>
        /// PopulationSizeTextbox_TextChanged handler
        /// </summary>
        /// <param name="sender"></param>
        /// <param name="e"></param>
        private void PopulationSizeTextBox_TextChanged(object sender, TextChangedEventArgs e)
        {
            try
            {
                Optimizer.PopulationSize = Int32.Parse(PopulationSizeTextBox.Text);
            }
            catch
            {
                ModernDialog.ShowMessage("Incorrenct value. Please, enter a number.", "Warning", MessageBoxButton.OK);
            }
        }


        /// <summary>
        /// ReproductionCyclesTextbox_TextChanged handler
        /// </summary>
        /// <param name="sender"></param>
        /// <param name="e"></param>
        private void ReproductionCyclesTextBox_TextChanged(object sender, TextChangedEventArgs e)
        {
            try
            {
                Optimizer.ReproductiveCyclesNumber = Int32.Parse(ReproductionCyclesTextBox.Text);
            }
            catch
            {
                ModernDialog.ShowMessage("Incorrenct value. Please, enter a number.", "Warning", MessageBoxButton.OK);
            }
        }

        /// <summary>
        /// TournamentSizeTextbox_TextChanged handler
        /// </summary>
        /// <param name="sender"></param>
        /// <param name="e"></param>
        private void TournamentSizeTextBox_TextChanged(object sender, TextChangedEventArgs e)
        {
            try
            {
                Optimizer.TournamentSize = Int32.Parse(TournamentSizeTextBox.Text);
            }
            catch
            {
                ModernDialog.ShowMessage("Incorrenct value. Please, enter a number.", "Warning", MessageBoxButton.OK);
            }
        }

        /// <summary>
        /// StopCriterionTextBox_TextChanged handler
        /// </summary>
        /// <param name="sender"></param>
        /// <param name="e"></param>
        private void StopCriterionTextBox_TextChanged(object sender, TextChangedEventArgs e)
        {
            try
            {
                Optimizer.StopCriterion = Int32.Parse(StopCriterionTextBox.Text);
            }
            catch
            {
                ModernDialog.ShowMessage("Incorrenct value. Please, enter a number.", "Warning", MessageBoxButton.OK);
            }
        }


        /// <summary>
        /// MutationProbabilityTextBox_TextChanged handler
        /// </summary>
        /// <param name="sender"></param>
        /// <param name="e"></param>
        private void MutationProbabilityTextBox_TextChanged(object sender, TextChangedEventArgs e)
        {
            try
            {
                Optimizer.MutationProbability = float.Parse(MutationProbabilityTextBox.Text, CultureInfo.InvariantCulture);
            }
            catch
            {
                ModernDialog.ShowMessage("Incorrenct value. Please, enter a number.", "Warning", MessageBoxButton.OK);
            }
        }

        /// <summary>
        /// CrossoverProbabilityTextBox_TextChanged handler
        /// </summary>
        /// <param name="sender"></param>
        /// <param name="e"></param>
        private void CrossoverProbabilityTextBox_TextChanged(object sender, TextChangedEventArgs e)
        {
            try
            {
                Optimizer.CrossoverProbability = float.Parse(CrossoverProbabilityTextBox.Text, CultureInfo.InvariantCulture);
            }
            catch
            {
                ModernDialog.ShowMessage("Incorrenct value. Please, enter a number.", "Warning", MessageBoxButton.OK);
            }
        }
        #endregion
    }
}
