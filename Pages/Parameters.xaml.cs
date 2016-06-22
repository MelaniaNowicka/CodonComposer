using CodonOptimizer.Classes;
using FirstFloor.ModernUI.Windows.Controls;
using Microsoft.VisualBasic.FileIO;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
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
            RestrictionSitesPoolListBox.ItemsSource = new List<string>();
            ImportEnzymes();
            RestrictionSitesPoolListBox.ItemsSource = SeqParser.enzymesToSequences.Keys.ToList();
            Optimizer.EnzymeSitesToRemoval = new List<string>();
            Optimizer.AHomopolymersRemoval = false;
            Optimizer.NcRestrictions = false;
            Optimizer.RestrEnzymeSitesToRemoval = false;
            RestrictionSitesPoolListBox.IsEnabled = false;
            RestrictionSitesToRemovalListBox.IsEnabled = false;
            AddToRemovalButton.IsEnabled = false;
            UndoAddToRemovalButton.IsEnabled = false;
        }

        /// <summary>
        /// Enzymes import from file.
        /// </summary>
        private void ImportEnzymes()
        {
            string path = System.IO.Path.Combine(Directory.GetCurrentDirectory(), "restren.csv");

            SeqParser.enzymesToSequences = new Dictionary<string, List<string>>();

            using (TextFieldParser parser = new TextFieldParser(path))
            {
                parser.SetDelimiters(new string[] { ";" });

                if (parser.EndOfData)
                {
                    // modern dialog initialization
                    string message = "The restriction enzymes file is empty.";
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

                        SeqParser.enzymesToSequences.Add(fields[0], tmpList);
                    }
                }
            }
        }

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

        /// <summary>
        /// AddToRemovalButton_Click handler
        /// </summary>
        /// <param name="sender"></param>
        /// <param name="e"></param>
        private void AddToRemovalButton_Click(object sender, RoutedEventArgs e)
        {
            foreach (var item in RestrictionSitesPoolListBox.SelectedItems)
            {
                if (!RestrictionSitesToRemovalListBox.Items.Contains(item))
                {
                    RestrictionSitesToRemovalListBox.Items.Add(item);
                    Optimizer.EnzymeSitesToRemoval.Add(item.ToString());
                }
            }

            foreach (string str in Optimizer.EnzymeSitesToRemoval)
            {
                Console.Write(str + " ");
            }
        }

        /// <summary>
        /// UndoAddToRemovalButton_Click handler
        /// </summary>
        /// <param name="sender"></param>
        /// <param name="e"></param>
        private void UndoAddToRemovalButton_Click(object sender, RoutedEventArgs e)
        {
            List<string> itemsToRemoval = new List<string>();

            foreach (string item in RestrictionSitesToRemovalListBox.SelectedItems)
            {
                itemsToRemoval.Add(item);
            }

            foreach (string item in itemsToRemoval)
            {
                RestrictionSitesToRemovalListBox.Items.Remove(item);
                Optimizer.EnzymeSitesToRemoval.Remove(item.ToString());
            }

            foreach (string str in Optimizer.EnzymeSitesToRemoval)
            {
                Console.Write(str + " ");
            }
        }
        #endregion

        private void AHomopolymersRemovalCheckBox_Checked(object sender, RoutedEventArgs e)
        {
            Optimizer.AHomopolymersRemoval = true;
        }

        private void AHomopolymersRemovalCheckBox_Unchecked(object sender, RoutedEventArgs e)
        {
            Optimizer.AHomopolymersRemoval = false;
        }

        private void NcOptimizationCheckBox_Checked(object sender, RoutedEventArgs e)
        {
            Optimizer.NcRestrictions = true;
        }

        private void NcOptimizationCheckBox_Unchecked(object sender, RoutedEventArgs e)
        {
            Optimizer.NcRestrictions = false;
        }

        private void RestrEnzymeSitesToRemovalCheckBox_Checked(object sender, RoutedEventArgs e)
        {
            RestrictionSitesPoolListBox.IsEnabled = true;
            RestrictionSitesToRemovalListBox.IsEnabled = true;
            AddToRemovalButton.IsEnabled = true;
            UndoAddToRemovalButton.IsEnabled = true;
            Optimizer.RestrEnzymeSitesToRemoval = true;
        }

        private void RestrEnzymeSitesToRemovalCheckBox_Unchecked(object sender, RoutedEventArgs e)
        {
            RestrictionSitesPoolListBox.IsEnabled = false;
            RestrictionSitesToRemovalListBox.IsEnabled = false;
            AddToRemovalButton.IsEnabled = false;
            UndoAddToRemovalButton.IsEnabled = false;
            Optimizer.RestrEnzymeSitesToRemoval = false;
        }

        private void NcRestrictionsCheckbox_Checked(object sender, RoutedEventArgs e)
        {
            Optimizer.NcRestrictions = true;
        }

        private void NcRestrictionsCheckbox_Unchecked(object sender, RoutedEventArgs e)
        {
            Optimizer.NcRestrictions = false;
        }

        private void MinimalNcTextBox_TextChanged(object sender, TextChangedEventArgs e)
        {
            try
            {
                Optimizer.MinimalNc = Convert.ToInt32(MinimalNcTextBox.Text);
            }
            catch
            {

            }
        }

        private void MaximalNcTextBox_TextChanged(object sender, TextChangedEventArgs e)
        {
            try
            {
                Optimizer.MaximalNc = Convert.ToInt32(MaximalNcTextBox.Text);
            }
            catch
            {

            }

        }


        private void NcRestrictionsMaintenanceCheckBox_Unchecked(object sender, RoutedEventArgs e)
        {
            Optimizer.MaintainOriginalNc = false;
        }

        private void NcRestrictionsMaintenanceCheckBox_Checked(object sender, RoutedEventArgs e)
        {
            Optimizer.MaintainOriginalNc = true;
            ModernDialog.ShowMessage("Nc maintenance require another version of algorithm. The results may be worse.", "Information", MessageBoxButton.OK);
        }
    }
}
