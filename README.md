# CRN Simulator

This repository enables quick simulation of a CRN (Chemical Reaction Network) that is specified using a JSON file.
It outputs the distribution of the different species through time.

It can interface with Stochkit (through the ```--stochkit_path``` flag) or CMEPy (using the ```--cmepy``` flag).

Example of usage:

```bash
python crn.py --model model.json --show_plots --cmepy
```

which should output:

![Screenshot](https://raw.githubusercontent.com/amandaprorok/crn/master/img/screenshot.png)

# Privacy

## A quantitative privacy model.

We are interested in securing the operation of robotic networks composed of heterogeneous agents. Since any given robot type plays a role that may be critical in guaranteeing continuous and failure-free operation of the system, it is beneficial to hide information that reveals the individual robot types and, thus, their roles. We propose a method that quantifies how easy it is for an adversary to identify the type of any of the robots, based on an outside observation of the system’s behavior. We draw from the theory of differential privacy, and develop an analytical model of the information leakage. This model allows us to analyze the privacy of the system, as its parameters vary.

## Relevant publications

This code was used to produce the results seen in the following privacy-related papers:

A. Prorok, V. Kumar, Towards Differentially Private Aggregation of Heterogeneous Robots, 13th International Symposium on Distributed Autonomous Robotic Systems (DARS), 2016. [PDF](http://prorok.me/wp-content/uploads/2015/02/DARS-2016_Prorok.pdf)

A. Prorok, V. Kumar, A Macroscopic Privacy Model for Heterogeneous Robot Swarms, International Conference on Swarm Intelligence, 2016, [PDF](http://prorok.me/wp-content/uploads/2016/06/2016_PrivateSwarms_ANTS.pdf)

A. Prorok, V. Kumar, A Macroscopic Model for Differential Privacy in Dynamic Robotic Networks. Swarm Intelligence, under review. [PDF (preprint)](http://prorok.me/wp-content/uploads/2016/12/SI_Prorok_R1.pdf)
