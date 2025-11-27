# Ttbar resonance
This repository contains the code needed to generate events for the ttbar/heavy-higgs interference.
This events are used to train an algorithm to regress the parameters that govern the BSM models.

## Event geneneration instructions

Please have a look [here](https://www.notion.so/MadGraph-generation-286570f8119b80c4a6a8d9cec861f13c?source=copy_link).

## HTCondor:

Run:
```
condor_submit -append "ARGS=Signal_c1p55 SIGNAL 1.5" Submit.sub
```