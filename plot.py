import os
import re
import argparse
import matplotlib.pyplot as plt

def plot_scores_by_iteration(input_dir, output_dir):
    """
    Reads the files in the specified directory, extracts score data,
    generates a separate plot for each score, and a combined plot for all scores,
    saving them in the specified output directory.

    Args:
        input_dir (str): Path to the directory containing the score files.
        output_dir (str): Directory to save the plots.
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Initialize a dictionary to collect score data and iterations
    score_data = {}
    iterations = []

    # Read each file in the directory
    for filename in os.listdir(input_dir):
        if filename.startswith("score_summary") and filename.endswith(".txt"):
            filepath = os.path.join(input_dir, filename)
            
            # Extract iteration number from the filename
            match = re.search(r"score_summary_system_\d+_(\d+)_traj_reimaged\.txt", filename)
            if match:
                iteration = int(match.group(1))
                iterations.append(iteration)

                with open(filepath, 'r') as f:
                    for line in f:
                        parts = line.strip().split()
                        for part in parts[1:]:
                            score_name, score_value = part.split(":")
                            score_value = float(score_value)
                            if score_name not in score_data:
                                score_data[score_name] = []
                            score_data[score_name].append((iteration, score_value))

    # Sort data by iteration and organize by each score
    for score_name in score_data:
        score_data[score_name] = sorted(score_data[score_name], key=lambda x: x[0])

    # Generate individual plots for each score
    for score_name, values in score_data.items():
        iterations, scores = zip(*values)
        
        # Create the plot for each score
        plt.figure(figsize=(8, 5))
        plt.plot(iterations, scores, marker='o', label=score_name)
        plt.xlabel("Iteration")
        plt.ylabel("Score Value")
        plt.title(f"Score Trend for {score_name.upper()} by Iteration")
        plt.grid(True)

        # Save the plot for the specific score
        output_path = os.path.join(output_dir, f"{score_name}_summary.png")
        plt.savefig(output_path)
        plt.close()

    # Generate a combined plot for all scores
    plt.figure(figsize=(10, 6))
    for score_name, values in score_data.items():
        iterations, scores = zip(*values)
        plt.plot(iterations, scores, marker='o', label=score_name)

    # Configure the combined plot
    plt.xlabel("Iteration")
    plt.ylabel("Score Value")
    plt.title("Combined Score Trend by Iteration")
    plt.grid(True)
    plt.legend()

    # Save the combined plot
    overall_output_path = os.path.join(output_dir, "overall_scores_summary.png")
    plt.savefig(overall_output_path)
    plt.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate score plots by iteration.")
    parser.add_argument("--input_dir", required=True, help="Path to the directory containing the score files.")
    parser.add_argument("--output_dir", required=True, help="Path to the directory where the plots will be saved.")
    args = parser.parse_args()

    plot_scores_by_iteration(args.input_dir, args.output_dir)
