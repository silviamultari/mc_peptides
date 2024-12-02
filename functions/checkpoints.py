from pathlib import Path
import pickle

def save_checkpoint(
    i, 
    checkpoint_file_base, 
    ref_seq, 
    old_pos, 
    old_i, 
    interval=10
):
    """
    Saves a checkpoint periodically based on the given interval.

    Args:
        i (int): The current iteration.
        checkpoint_file_base (Path): Base directory or file path for checkpoints.
        ref_seq (object): Reference sequence to save.
        old_pos (object): Old position data to save.
        old_i (object): Previous index data to save.
        interval (int): Interval at which to save checkpoints (default: 10).
    """
    if i % interval == 0:
        # Generate the checkpoint filename
        checkpoint_filename = checkpoint_file_base / f"checkpoint_{i}.pkl"
        
        # Prepare the checkpoint data
        checkpoint_data = {
            'ref_seq': ref_seq,
            'old_pos': old_pos,
            'old_i': old_i,
            'iteration': i + 1
        }
        
        # Save the checkpoint using pickle
        with open(checkpoint_filename, 'wb') as file:
            pickle.dump(checkpoint_data, file)
        
        print(f"Checkpoint saved: {checkpoint_filename}")
