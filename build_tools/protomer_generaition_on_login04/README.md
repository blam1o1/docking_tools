## Protonation on login04

Current implementation of protonation workflow is unable to restrain job to a single thread. 
Therefore you are not able to run protonation jobs in parallel on a single node. This submission 
strategy automatically submits protomer generation jobs across all available nodes making sure not to 
submit to the same node twice. Scheduler on each node will loop through batches of 50K smiles per job.
