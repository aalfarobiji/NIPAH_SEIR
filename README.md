# NIPAH_SEIR
Mathematical modeling of Nipah virus transmission using SEIR-D under a One Health framework. The model compares Bat–Pig–Human and Bat–Human spillover pathways, evaluates outbreak thresholds (R₀), and explores intervention scenarios to support outbreak preparedness and risk assessment in Indonesia.

Nipah Virus Transmission Modeling in Indonesia
This repository presents a mathematical modeling framework to analyze potential Nipah virus (NiV) outbreak scenarios in Indonesia using a SEIR-D compartmental approach within a One Health framework. The model explores two major zoonotic transmission pathways and evaluates outbreak dynamics, key epidemiological parameters, and intervention strategies.

Background
Nipah virus is a highly lethal zoonotic pathogen with a case fatality rate (CFR) of 40–75%, recognized by the WHO as a potential pandemic threat. The natural reservoir of the virus is fruit bats (Pteropus spp.), which are widely distributed across Indonesia. Spillover transmission may occur either directly from bats to humans or through intermediate hosts such as pigs.

Modeling Approach
The study uses deterministic compartmental models (SEIR-D) to simulate disease dynamics across multiple host populations. Two transmission pathways are considered:

1. Bat–Pig–Human (BPH) Pathway
* Three interacting populations: bats (reservoir), pigs (amplifier host), and humans (final host).
* Five transmission routes:
bat-to-bat
bat-to-pig
pig-to-pig
pig-to-human
human-to-human

2. Bat–Human (BH) Pathway
* Two populations: bats and humans.
* Transmission routes:
bat-to-bat
bat-to-human spillover
human-to-human transmission.

Both models track transitions between Susceptible, Exposed, Infectious, Recovered, and Dead (SEIR-D) compartments.
Key Parameter
A critical parameter in the model is μ (mu), representing the probability that exposed individuals progress to infectious clinical cases. Simulation scenarios include:

μ value	    Estimated R₀	    Interpretation
0.1	        ~0.42	            Limited transmission; outbreak dies out
0.3        	~1.26	            Sustained transmission possible
0.5	        ~2.10	            Rapid epidemic expansion

The epidemic threshold occurs at μ ≈ 0.24 (R₀ = 1). 

Key Findings

Amplifier role of pigs:
1. In the Bat–Pig–Human pathway, pigs significantly amplify infection dynamics before human cases emerge.
2. Sensitivity to μ:
Increasing μ accelerates outbreak growth and increases cumulative mortality.
3. High mortality burden:
With CFR ≈ 71%, even moderate outbreaks can produce severe public health impacts.
4. Intervention potential:
Combined interventions reducing human transmission (β_HH) and infection progression (μ) can lower R₀ from ~1.26 to ~0.91, preventing sustained outbreaks.

Policy Implications (One Health)
The results highlight the importance of integrated surveillance and preparedness across humans, livestock, and wildlife, including:
* biosurveillance of Pteropus bat populations
* biosecurity in pig farms in high-risk areas
* early detection of atypical pneumonia or encephalitis cases
* rapid isolation and contact tracing
* strengthened laboratory capacity for Nipah PCR diagnostics.

Purpose of the Repository
This repository aims to:
* provide a modeling framework for zoonotic outbreak analysis
* support risk assessment and preparedness planning
* explore scenario simulations for Nipah virus emergence in Indonesia
* contribute to evidence-based One Health policy development
