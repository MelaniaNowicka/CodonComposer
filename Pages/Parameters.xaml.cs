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
            RestrictionSitesPoolListBox.ItemsSource = new List<string>();
            RestrictionSitesPoolListBox.ItemsSource = SeqParser.enzymesToSequences.Keys.ToList();
            Optimizer.EnzymeSitesToRemoval = new List<string>();
            Optimizer.AHomopolymersRemoval = false;
            Optimizer.NTerminusOptimization = false;
            Optimizer.RestrEnzymeSitesToRemoval = false;
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

        private void NTerminusOptimizationCheckBox_Checked(object sender, RoutedEventArgs e)
        {
            Optimizer.NTerminusOptimization = true;
        }

        private void NTerminusOptimizationCheckBox_Unchecked(object sender, RoutedEventArgs e)
        {
            Optimizer.NTerminusOptimization = false;
        }

        private void RestrEnzymeSitesToRemovalCheckBox_Checked(object sender, RoutedEventArgs e)
        {
            Optimizer.RestrEnzymeSitesToRemoval = true;
        }

        private void RestrEnzymeSitesToRemovalCheckBox_Unchecked(object sender, RoutedEventArgs e)
        {
            Optimizer.RestrEnzymeSitesToRemoval = false;
        }


    }
}
