"""
Dreemium Messages of Success

Print messages to indicate successful completion of one of Dreemium's programs.
"""

import random


def dms(func):
    def wrapper(*args, **kwargs):
        ret = func(*args, **kwargs)
        print(MessagePicker().quote)
        return ret
    return wrapper


class MessagePicker(object):
    def __init__(self):
        self.num = random.randint(1, len(self._messages) - 2)

    @property
    def message(self):
        return self._messages[self.num].strip()

    @property
    def quote(self):
        return f"Dreemium Message of Success #{self.num}: '{self.message}'"

    _messages = """
    Detect Mutations Speedily
    Deconvolute Multiple Structures
    Don't MaP with SHAPE
    Dancing at Midnight to Salsa
    Asp-Met-Ser
    Dreams, Marianne (by Storr)
    Doing Minipreps on a Saturday
    Debugging My Software
    Dreaming of a Manuscript in Science
    Drosophila Melanogaster Sevenless
    Doping Makes Semiconductors
    Deprived of Much Sleep
    Drop Mentos into Soda
    Department of Materials Science
    Dijon Mustard Smacks
    DMSO: Moistureless Solvent
    Destroy Microbes with Sulfonamides
    Defected MIT Student
    Drug Many Stem-loops
    Duststorms Move Soil
    Dark Monolith in Space
    Daytime Megawatts from Solar
    Despised Murder-hornets Sting
    Delusional Misidentification Syndrome
    Dry-ice Melts not: Sublimates
    Diatoms Make Shells
    Driving in Massachusetts Sucks
    Dissertation-Making is Stressful
    Deep Marine Sea-vents
    Dielectric Materials Shield
    Do Meticulous Science
    Deceit Makes Secrets
    DNA Multiplexed Sequencing
    Dating Millions of Samples
    Deduction Master Sherlock
    Document Management System
    Degrees° Minutes' Seconds"
    Details Matter Significantly
    Diligence is Mandatory for Success
    Decapitated Moving-trucks got Storrowed
    Dreading a Major Scoop
    Dextrose: Munch Sweetly
    Dying of Med-student Syndrome
    Diagnostic Medical Sonography
    Delivering Mail by Snail
    Dairy MilkShake
    Delaware Makes Startups
    Devoured by Mackerel Sharks
    Dumping Mountains of Snow
    Don't Miss the Sunset
    Drenching Monsoon Season
    David: a Michelangelo Sculpture
    Delicious Midnight Snacking
    Distrust Most Structures
    Drowsily, Melatonin Slumbers
    Driver Monitoring System
    Disco Music: will Survive
    Driving on Martian Soil
    Doubting My Simulations
    Dashi Miso Soup
    Detect Methionine with Sulfur-35
    Division of Medical Sciences
    Dinner of Meatballs and Spaghetti
    Dear Mansplainer: Silence!
    Degree: Master of Science
    Darn Mangled Sweaters
    A/G/U A/C C/G
    Determine MHV Structure
    Don Masks Safely
    🎼 Do-Mi-So 🎶
    Dynamic Message Sign
    Dinosaurs: Mesozoic Sauropods
    Daily Mismatched Socks
    Disrupt Membranes with SDS
    Direct Message me on Slack
    Data from Money is all of Science
    Drive Me to Singapore
    Disrupting Markets with Software
    Double My Stipend (please)
    Duly, Mimes are Silent
    Do Mitochondria Swim?
    Darwin Memorized Species
    Deus ex Machina in Sophocles
    Duckling Matures into Swan
    Downstairs, Monsters Slumber
    Distribute Modules with Setuptools
    Diogenes Mocked Society
    Defuse Mines Safely
    Davy's Minedamp Safety-lamp
    Don't Mash Spiders
    Dragon Marauder Smaug
    Doge Meme $
    Dodos from Mare aux Songes
    Dido Married Sychaeus
    Diseased Martians Succumbed
    Don't Monger Stereotypes
    DaVinci Made Sketches
    Darth Maul Survived
    Decorated Michelangelo, the Sistine
    Duet of Matter and Spacetime
    DREEM Me up, Scotty!
    Diagnose Methane using Sulfhydryls
    Darkened Marshmallows Smolder
    Develop Machine-learning Solutions
    Daily Meditation for the Soul
    Desmond's Melodious Saxophone
    Daring Motorcycle Stunts
    Demeter Moves the Seasons
    Desert in Mexico: Sonora
    DuPont Manufacturers Sorona
    Diphenhydramine Makes you Sleepy
    Decompose Matrices Spectrally
    Deeds Manifest Sentiments
    Dawn -> Midday -> Sunset
    Demise of Mozart: Streptococci?
    D Minor Sonata
    Duplicating Manuscripts via Selenium
    Deep Mariana Submersible
    Doomscrolling Magnifies Stress
    Digital Microprocessors: Sentient?
    Doomsday Machine in Strangelove
    Dead Man's Switch
    """.split("\n")
