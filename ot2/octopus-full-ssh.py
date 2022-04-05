import sys
from opentrons import protocol_api

metadata = {
    'apiLevel': '2.12',
    'protocolName': 'OCTOPUS OT-2 Protocol: SSH, All Stages, Offset Calibrated',
    'author': 'Octant',
    'description': 'OT-2 Protocol for running Twist 96-Plex Library Kit (formally iGenomX RIPTIDE) preparation adapted for OCTOPUS. This version is for running on the OT-2 directly through SSH, including offset calibration, and runs through all stages.'
}

def run(protocol: protocol_api.ProtocolContext):
    ### specify the number of plates to process: integer between 1 and 8 inclusive
    number_of_plates = 0
    ### specify True if running protocol on OT-2 command line, otherwise False
    ssh_mode = True

    protocol.set_rail_lights(True)
    Run = OCTOPUS(protocol, ssh_mode)
    Run.Pool(number_of_plates)
    Run.Post_Pool()
    protocol.set_rail_lights(False)

class OCTOPUS:
    def __init__(self, protocol, ssh_mode=False):
        ### tip-rack definitions
            # opentrons 200 ul filtered tips
        tips_200_1 = protocol.load_labware('opentrons_96_filtertiprack_200ul', '7')
        tips_200_2 = protocol.load_labware('opentrons_96_filtertiprack_200ul', '8')
        tips_200_3 = protocol.load_labware('opentrons_96_filtertiprack_200ul', '10')
        tips_200_4 = protocol.load_labware('opentrons_96_filtertiprack_200ul', '11')
            # opentrons 20 ul filtered tips
        tips_20_1 = protocol.load_labware('opentrons_96_filtertiprack_20ul', '3') 

        if ssh_mode:
            tips_200_1.set_offset(x=-0.60, y=1.00, z=-0.20)
            tips_200_2.set_offset(x=-0.60, y=1.00, z=-0.20)
            tips_200_3.set_offset(x=-0.60, y=1.00, z=-0.20)
            tips_200_4.set_offset(x=-0.60, y=1.00, z=-0.20)
            tips_20_1.set_offset(x=0.00, y=1.00, z=-0.10)

        self.tips_200 = [tips_200_1, tips_200_2, tips_200_3, tips_200_4]
        self.tips_20 = [tips_20_1]

        ### pipette definitions
        self.p300 = protocol.load_instrument('p300_multi_gen2', 'left', tip_racks=self.tips_200)
        self.p20 = protocol.load_instrument('p20_multi_gen2', 'right', tip_racks=self.tips_20)

        ### load magnetic module
        self.mag_mod = protocol.load_module('magnetic module gen2', '9')

        ### pass down protocol api throughout object
        self.protocol = protocol

        ### pass down ssh-mode status throughout object
        self.ssh_mode = ssh_mode

    def Pool(self, num_plates=0):
        NUM_PLATES = num_plates
        if NUM_PLATES < 1 or NUM_PLATES > 8:
            raise Exception("Error: Number of plates must be between 1 and 8 inclusive.")

        ### labware specific to Pool stage
        self.pool_plate = self.mag_mod.load_labware('nest_96_wellplate_2ml_deep', label='Pool Plate')
        self.dna_plate = self.protocol.load_labware('biorad_96_wellplate_200ul_pcr', '6', label='DNA Plate')
        self.consolidation_plate = self.protocol.load_labware('biorad_96_wellplate_200ul_pcr', '4', label='Consolidation Plate')
        self.pre_pool_plate_1 = self.protocol.load_labware('appliedbiosystems_384_wellplate_40ul', '1', label='Pre-pool Plate 1 (1-4)')
        self.pre_pool_plate_2 = self.protocol.load_labware('appliedbiosystems_384_wellplate_40ul', '2', label='Pre-pool Plate 2 (5-8)')

        SSH_MODE = self.ssh_mode
        if SSH_MODE:
            self.pool_plate.set_offset(x=0.10, y=1.60, z=0.00)
            self.dna_plate.set_offset(x=0.40, y=1.10, z=0.10)
            self.consolidation_plate.set_offset(x=0.30, y=1.10, z=0.00)
            self.pre_pool_plate_1.set_offset(x=0.00, y=0.80, z=0.60)
            self.pre_pool_plate_2.set_offset(x=0.00, y=0.80, z=0.60)

        ### reagent locations
        DNA_pool = self.pool_plate.wells_by_name()['A1']
        Purification_Beads = self.pool_plate.wells_by_name()['A3']
        Purification_Beads_discard = self.pool_plate.wells_by_name()['A4']
        Ethanol = [self.pool_plate.wells_by_name()['A6'], self.pool_plate.wells_by_name()['A9']]
        Ethanol_discard = [self.pool_plate.wells_by_name()['A7'], self.pool_plate.wells_by_name()['A10']]
        Tris_HCl = self.pool_plate.wells_by_name()['A12']

        ### consolidate Reaction A plate after iGenomX_1 into consolidation plate
        self.set_flow_rate(48, 48, 48)
            # consolidating VE plates 1-4 to Pool Plate columns 1-4
        if NUM_PLATES >= 1:
            self.p300.consolidate(5, [self.pre_pool_plate_1.rows_by_name()['A'][num_plate].bottom(z=0.5) for num_plate in range(0, 23, 2)], self.consolidation_plate.columns('1'))
        if NUM_PLATES >= 2:
            self.p300.consolidate(5, [self.pre_pool_plate_1.rows_by_name()['A'][num_plate].bottom(z=0.5) for num_plate in range(1, 24, 2)], self.consolidation_plate.columns('2'))
        if NUM_PLATES >= 3:
            self.p300.consolidate(5, [self.pre_pool_plate_1.rows_by_name()['B'][num_plate].bottom(z=0.5) for num_plate in range(0, 23, 2)], self.consolidation_plate.columns('3'))
        if NUM_PLATES >= 4:
            self.p300.consolidate(5, [self.pre_pool_plate_1.rows_by_name()['B'][num_plate].bottom(z=0.5) for num_plate in range(1, 24, 2)], self.consolidation_plate.columns('4'))
            # consolidating VE plates 5-8 to Pool Plate columns 5-8
        if NUM_PLATES >= 5:
            self.p300.consolidate(5, [self.pre_pool_plate_2.rows_by_name()['A'][num_plate].bottom(z=0.5) for num_plate in range(0, 23, 2)], self.consolidation_plate.columns('5'))
        if NUM_PLATES >= 6:
            self.p300.consolidate(5, [self.pre_pool_plate_2.rows_by_name()['A'][num_plate].bottom(z=0.5) for num_plate in range(1, 24, 2)], self.consolidation_plate.columns('6'))
        if NUM_PLATES >= 7:
            self.p300.consolidate(5, [self.pre_pool_plate_2.rows_by_name()['B'][num_plate].bottom(z=0.5) for num_plate in range(0, 23, 2)], self.consolidation_plate.columns('7'))
        if NUM_PLATES == 8:
            self.p300.consolidate(5, [self.pre_pool_plate_2.rows_by_name()['B'][num_plate].bottom(z=0.5) for num_plate in range(1, 24, 2)], self.consolidation_plate.columns('8'))

        ### further consolidate plate columns into single wells on the DNA pool deep well plate
        self.set_flow_rate(48*2, 48, 48)
        for i in range(NUM_PLATES):
            self.p300.pick_up_tip(self.tips_200[0].rows()[7-i][NUM_PLATES])
            # aspirate first four wells
            for j in range(4):
                self.p300.aspirate(50, self.consolidation_plate.rows()[j][i].bottom(z=0.5))
            self.p300.dispense(200, self.pool_plate.rows()[i][0])
            # aspirate last four wells
            for j in range(4):
                self.p300.aspirate(50, self.consolidation_plate.rows()[j+4][i].bottom(z=0.5))
            self.p300.dispense(200, self.pool_plate.rows()[i][0])
            self.p300.drop_tip()

        ### add and incubate Purification Beads with pooled samples
        self.p300.pick_up_tip()
        self.set_flow_rate(48*4, 48*4, 48)
        self.p300.mix(10, 200, Purification_Beads.bottom(z=10.0))
            # adding beads to pool
        self.set_flow_rate(48*2, 48*4, 48)
        for i in range(3):
            self.p300.aspirate(200, Purification_Beads)
            self.p300.dispense(200, DNA_pool)
        self.p300.aspirate(120, Purification_Beads)
        self.p300.dispense(120, DNA_pool)
            # mix beads well with pool and incubate
        self.set_flow_rate(48*4, 48*6, 48)
        self.p300.mix(10, 200, DNA_pool.bottom(z=10.0))
        self.set_flow_rate(48*4, 48*4, 48)
        for i in range(20):
            self.p300.aspirate(200, DNA_pool.bottom(z=10.0))
            self.p300.dispense(200, DNA_pool.center())
        self.set_flow_rate(48*4, 48*6, 48/6)
        self.p300.mix(20, 200, DNA_pool.bottom(z=10.0))
        self.p300.move_to(DNA_pool.top())
        self.protocol.delay(seconds=5.0)
        self.p300.blow_out(DNA_pool.top())
        self.p300.drop_tip()
        self.protocol.delay(minutes=10)
        self.mag_mod.engage()
        self.protocol.delay(minutes=12)

        ### discard supernatant in stages to minimize bead loss
        self.set_flow_rate(48/2, 48, 48)
        self.p300.pick_up_tip()
        self.p300.aspirate(200, DNA_pool.bottom(z=0.5))
        self.p300.dispense(200, Purification_Beads_discard)
        self.p300.move_to(Purification_Beads_discard.top(z=10.0))
        self.protocol.delay(minutes=4)
        self.p300.aspirate(200, DNA_pool.bottom(z=0.5))
        self.p300.dispense(200, Purification_Beads_discard)
        self.p300.move_to(Purification_Beads_discard.top(z=10.0))
        self.protocol.delay(minutes=2)
        self.p300.aspirate(200, DNA_pool.bottom(z=0.5))
        self.p300.dispense(200, Purification_Beads_discard)
        self.p300.move_to(Purification_Beads_discard.top(z=10.0))
        self.protocol.delay(minutes=2)
        self.p300.aspirate(200, DNA_pool.bottom(z=0.5))
        self.p300.dispense(200, Purification_Beads_discard)
        self.p300.move_to(Purification_Beads_discard.top(z=10.0))
        self.protocol.delay(minutes=2)
        self.p300.aspirate(200, DNA_pool.bottom(z=0.5))
        self.p300.dispense(200, Purification_Beads_discard)
        self.p300.move_to(Purification_Beads_discard.top(z=10.0))
        self.protocol.delay(minutes=2)
        self.p300.aspirate(120, DNA_pool.bottom(z=0.5))
        self.p300.dispense(120, Purification_Beads_discard)
        self.p300.drop_tip()
        self.set_flow_rate(5, 48, 48)
        self.p20.pick_up_tip()
        self.p20.aspirate(20, DNA_pool.bottom())
        self.p20.drop_tip()

        ### wash with 80% ethanol twice
        self.set_flow_rate(48, 48, 48)
        self.p300.pick_up_tip()
        for i in range(6):
            self.p300.aspirate(200, Ethanol[0])
            self.p300.dispense(200, DNA_pool.top(z=-10.0))
        self.p300.move_to(DNA_pool.top(z=10.0))
        self.protocol.delay(seconds=40)
        self.discard_supernatant(1200, DNA_pool, Ethanol_discard[0].top(z=-10.0))
            # second time
        self.set_flow_rate(48, 48, 48)
        self.p300.pick_up_tip()
        for i in range(6):
            self.p300.aspirate(200, Ethanol[1])
            self.p300.dispense(200, DNA_pool.top(z=-10.0))
        self.p300.move_to(DNA_pool.top(z=10.0))
        self.protocol.delay(seconds=40)
        self.discard_supernatant(1200, DNA_pool, Ethanol_discard[1].top(z=-10.0))

        ### discard residual ethanol and let dry in stages to maximize ethanol removal
        self.set_flow_rate(5, 48, 48)
        self.p20.pick_up_tip()
        self.p20.aspirate(20, DNA_pool.bottom(z=-0.5))
        self.p20.drop_tip()
        self.protocol.delay(minutes=10)
        self.set_flow_rate(1, 48, 48)
        self.p20.pick_up_tip()
        self.p20.aspirate(5, DNA_pool.bottom(z=-0.5))
        self.p20.drop_tip()
        self.protocol.delay(minutes=5)
        self.set_flow_rate(1, 48, 48)
        self.p20.pick_up_tip()
        self.p20.aspirate(5, DNA_pool.bottom(z=-1.0))
        self.p20.drop_tip()

        ### resuspend beads well in 10 mM Tris-HCl pH 8
        self.mag_mod.disengage()
        self.set_flow_rate(48, 48, 48)
        self.p300.pick_up_tip()
        self.p300.aspirate(50, Tris_HCl.bottom(z=0.5))
        self.p300.dispense(50, DNA_pool)
        self.set_flow_rate(48, 48*2, 48/2)
        self.p300.mix(20, 30, DNA_pool.bottom(z=1.0))
        self.set_flow_rate(48/4, 48/4, 48/2)
        self.p300.mix(5, 45, DNA_pool.bottom(z=1.0))
        self.set_flow_rate(48, 48, 48)
        for i in range(20):
            self.p300.aspirate(30, DNA_pool.bottom(z=1.0))
            self.p300.dispense(30, DNA_pool.bottom(z=5.0))
        self.set_flow_rate(48/4, 48/4, 48/2)
        self.p300.mix(5, 45, DNA_pool.bottom(z=1.0))
        self.set_flow_rate(48, 48*2, 48/2)
        self.p300.mix(30, 30, DNA_pool.bottom(z=1.0))
        self.p300.blow_out(DNA_pool.bottom(z=5.0))
        self.p300.move_to(DNA_pool.top(z=10.0))
        self.protocol.delay(minutes=5)
            # resuspend settled beads
        self.set_flow_rate(48/4, 48/4, 48/2)
        self.p300.mix(5, 45, DNA_pool.bottom(z=1.0))
        self.set_flow_rate(48, 48, 48)
        for i in range(10):
            self.p300.aspirate(30, DNA_pool.bottom(z=1.0))
            self.p300.dispense(30, DNA_pool.bottom(z=5.0))
        self.set_flow_rate(48/4, 48/4, 48/2)
        self.p300.mix(5, 45, DNA_pool.bottom(z=1.0))
        self.set_flow_rate(48, 48*2, 48/2)
        self.p300.mix(30, 30, DNA_pool.bottom(z=1.0))
        self.p300.blow_out(DNA_pool.bottom(z=5.0))
        self.p300.drop_tip()
        self.protocol.delay(minutes=5)
        self.mag_mod.engage()
        self.protocol.delay(minutes=4)

        ### end of Pool Stage, prepare for Post-pool Stage
        self.set_flow_rate(48/4, 48, 48)
        self.p300.transfer(50, DNA_pool, self.dna_plate.wells_by_name()['A1'], blow_out=True, blowout_location='destination well')
        self.set_flow_rate(1, 5, 5)
        self.p20.pick_up_tip()
        self.p20.aspirate(10, DNA_pool.bottom())
        self.p20.dispense(5, self.dna_plate.wells_by_name()['A1'].top(z=-5.0))
        self.p20.dispense(5, self.dna_plate.wells_by_name()['A1'])
        self.p20.blow_out(self.dna_plate.wells_by_name()['A1'].top())
        self.p20.drop_tip()
        if SSH_MODE:
            input('Execute thermocycler protocol "iGenomX_Pre" on the DNA Plate. Replace former DNA Plate slot with plate containing Post-pool reagents. Empty trash contents if necessary. Discard deep-well plate. Afterward, spin down and secure DNA Plate on the Magnetic Module. Resume protocol for Post-pool stage. (Press ENTER)')
            input("ARE YOU SURE? Press ENTER to continue.")
        else:
            self.protocol.pause('Execute thermocycler protocol "iGenomX_Pre" on the DNA Plate. Replace former DNA Plate slot with plate containing Post-pool reagents. Empty trash contents if necessary. Discard deep-well plate. Afterward, spin down and secure DNA Plate on the Magnetic Module. Resume protocol for Post-pool stage.')
        self.mag_mod.disengage()
        self.protocol.home()

    def Post_Pool(self):
        ### redefine labware placements specific to Post-pool stage
        self.protocol.deck['9']._labware = None
        self.dna_plate = self.mag_mod.load_labware('biorad_96_wellplate_200ul_pcr', label='DNA Plate')
        del self.protocol.deck['6']
        self.post_pool_reagents = self.protocol.load_labware('biorad_96_wellplate_200ul_pcr', '6', label='Post-pool Reagents Plate')

        SSH_MODE = self.ssh_mode
        if SSH_MODE:
            self.dna_plate.set_offset(x=0.40, y=1.20, z=0.20)
            self.post_pool_reagents.set_offset(x=0.60, y=0.90, z=0.10)

        ### partition Post-pool stage into substages
        self.DNA_Capture_Library_Conversion()
        self.Library_Amplification()
        self.Size_Selection()

    def DNA_Capture_Library_Conversion(self):
        SSH_MODE = self.ssh_mode

        ### reagent locations
        Capture_Beads = self.post_pool_reagents.wells_by_name()['A1']
        HS_Buffer = self.post_pool_reagents.wells_by_name()['A2']
        NaOH = self.post_pool_reagents.wells_by_name()['A3']
        Bead_Wash_Buffer = [self.post_pool_reagents.wells_by_name()['A4'], self.post_pool_reagents.wells_by_name()['A4'], self.post_pool_reagents.wells_by_name()['A5']]

        ### prepare Capture Beads
        if SSH_MODE:
            input('Manually add 20 ul Capture Beads to column 2 of DNA Plate')
            input("ARE YOU SURE? Press ENTER to continue.")
        else:
            self.protocol.pause('Manually add 20 ul Capture Beads to column 2 of DNA Plate')
        self.mag_mod.engage()
        self.protocol.delay(minutes=4)
        #self.discard_supernatant(20, self.dna_plate.wells_by_name()['A2'])
        self.set_flow_rate(48/6, 48, 48)
        self.p300.pick_up_tip()
        self.p300.aspirate(20, self.dna_plate.wells_by_name()['A2'].bottom(z=0.5))
        self.p300.drop_tip()
        self.mag_mod.disengage()
            # wash with HS buffer
        self.set_flow_rate(48*2, 48*4, 48)
        self.p300.transfer(100, HS_Buffer, self.dna_plate.wells_by_name()['A2'], blow_out=True, blowout_location='destination well', mix_after=(10,100))
        self.mag_mod.engage()
        self.protocol.delay(minutes=10)
        self.discard_supernatant(100, self.dna_plate.wells_by_name()['A2'])
        self.mag_mod.disengage()
            # resuspend in HS buffer
        self.set_flow_rate(48*2, 48*4, 48)
        self.p300.transfer(20, HS_Buffer, self.dna_plate.wells_by_name()['A2'], blow_out=True, blowout_location='destination well', mix_after=(10,15))

        ### add pool samples into the Capture Beads
        self.set_flow_rate(48/4, 48, 48)
        self.p300.pick_up_tip()
        self.p300.aspirate(50, self.dna_plate.wells_by_name()['A1'].bottom(z=0.5))
        self.p300.dispense(50, self.dna_plate.wells_by_name()['A2'])
        self.set_flow_rate(48*2, 48*4, 48)
        self.p300.mix(10, 65, self.dna_plate.wells_by_name()['A2'])
        self.p300.blow_out(self.dna_plate.wells_by_name()['A2'].top())
        self.p300.move_to(self.dna_plate.wells_by_name()['A2'].top(z=10.0))
        self.protocol.delay(minutes=10)
            # resuspend same well again
        self.set_flow_rate(48*2, 48*4, 48)
        self.p300.mix(10, 65, self.dna_plate.wells_by_name()['A2'])
        self.p300.blow_out(self.dna_plate.wells_by_name()['A2'].top())
        self.p300.drop_tip()
        self.protocol.delay(minutes=10)
        self.mag_mod.engage()
        self.protocol.delay(minutes=10)
        self.discard_supernatant(70, self.dna_plate.wells_by_name()['A2'])
        self.mag_mod.disengage()

        ### denature with NaOH
        self.set_flow_rate(48*2, 48*4, 48)
        self.p300.transfer(100, NaOH, self.dna_plate.wells_by_name()['A2'], blow_out=True, blowout_location='destination well', mix_before=(3,95), mix_after=(10,95))
        self.protocol.delay(minutes=4)
        self.mag_mod.engage()
        self.protocol.delay(minutes=5)
        self.discard_supernatant(100, self.dna_plate.wells_by_name()['A2'])
        self.mag_mod.disengage()
        self.set_flow_rate(48*2, 48*4, 48)

        ### wash with Bead Wash Buffer three times
        self.p300.transfer(100, Bead_Wash_Buffer[0], self.dna_plate.wells_by_name()['A2'], blow_out=True, blowout_location='destination well', mix_after=(10,95))
        self.mag_mod.engage()
        self.protocol.delay(minutes=4)
        self.discard_supernatant(100, self.dna_plate.wells_by_name()['A2'])
        self.mag_mod.disengage()
            # second time
        self.set_flow_rate(48*2, 48*4, 48)
        self.p300.transfer(100, Bead_Wash_Buffer[1], self.dna_plate.wells_by_name()['A2'], blow_out=True, blowout_location='destination well', mix_after=(10,95))
        self.mag_mod.engage()
        self.protocol.delay(minutes=4)
        self.discard_supernatant(100, self.dna_plate.wells_by_name()['A2'])
        self.mag_mod.disengage()
            # third time
        self.set_flow_rate(48*2, 48*4, 48)
        self.p300.transfer(100, Bead_Wash_Buffer[2], self.dna_plate.wells_by_name()['A2'], blow_out=True, blowout_location='destination well', mix_after=(10,95))
        self.mag_mod.engage(height=11.0)
        self.protocol.delay(minutes=4)
        if SSH_MODE:
            input('Resume to complete wash step.')
            input("ARE YOU SURE? Press ENTER to continue.")
        else:
            self.protocol.pause('Resume to complete wash step.')
        self.discard_supernatant(100, self.dna_plate.wells_by_name()['A2'])
        self.discard_residual(self.dna_plate.wells_by_name()['A2'])

        ### proceed with Reaction B 
        if SSH_MODE:
            input('Resuspend each well in 18.5 ul Reaction B mix (without enzyme). Then, mix 0.5 ul Enzyme II into each well. Execute thermocycler protocol "iGenomX_2". Afterward, secure DNA Plate back on the Magnetic Module. Resume protocol. (Press ENTER)')
            input("ARE YOU SURE? Press ENTER to continue.")
        else:
            self.protocol.pause('Resuspend each well in 18.5 ul Reaction B mix (without enzyme). Then, mix 0.5 ul Enzyme II into each well. Execute thermocycler protocol "iGenomX_2". Afterward, secure DNA Plate back on the Magnetic Module. Resume protocol.')
        self.mag_mod.disengage()
        self.protocol.home()

    def Library_Amplification(self):
        SSH_MODE = self.ssh_mode

        ### reagent locations
        Bead_Wash_Buffer = [self.post_pool_reagents.wells_by_name()['A6'], self.post_pool_reagents.wells_by_name()['A6'], self.post_pool_reagents.wells_by_name()['A7']]

        ### resuspend beads after coming off thermocycler and discard reaction
        self.set_flow_rate(48*2, 48*6, 48)
        self.p300.pick_up_tip()
        self.p300.mix(10, 15, self.dna_plate.wells_by_name()['A2'])
        self.p300.blow_out(self.dna_plate.wells_by_name()['A2'].top())
        self.p300.drop_tip()
        self.mag_mod.engage()
        self.protocol.delay(minutes=4)
        self.discard_supernatant(19, self.dna_plate.wells_by_name()['A2'])
        self.mag_mod.disengage()

        ### wash with Bead Wash Buffer three times
        self.set_flow_rate(48*2, 48*4, 48)
        self.p300.transfer(100, Bead_Wash_Buffer[0], self.dna_plate.wells_by_name()['A2'], blow_out=True, blowout_location='destination well', mix_after=(10,95))
        self.mag_mod.engage()
        self.protocol.delay(minutes=4)
        self.discard_supernatant(100, self.dna_plate.wells_by_name()['A2'])
        self.mag_mod.disengage()
            # second time
        self.set_flow_rate(48*2, 48*4, 48)
        self.p300.transfer(100, Bead_Wash_Buffer[1], self.dna_plate.wells_by_name()['A2'], blow_out=True, blowout_location='destination well', mix_after=(10,95))
        self.mag_mod.engage()
        self.protocol.delay(minutes=4)
        self.discard_supernatant(100, self.dna_plate.wells_by_name()['A2'])
        self.mag_mod.disengage()
            # third time
        self.set_flow_rate(48*2, 48*4, 48)
        self.p300.transfer(100, Bead_Wash_Buffer[2], self.dna_plate.wells_by_name()['A2'], blow_out=True, blowout_location='destination well', mix_after=(10,95))
        self.mag_mod.engage(height=11.0)
        self.protocol.delay(minutes=4)
        if SSH_MODE:
            input('Resume to complete wash step.')
            input("ARE YOU SURE? Press ENTER to continue.")
        else:
            self.protocol.pause('Resume to complete wash step.')
        self.discard_supernatant(100, self.dna_plate.wells_by_name()['A2'])
        self.discard_residual(self.dna_plate.wells_by_name()['A2'])

        ### proceed with Amplification "Reaction C"
        if SSH_MODE:
            input('Resuspend each well in 50 ul Amplification "Reaction C" mix (primers added). Be sure to record which wells received which index primer. Execute thermocycler protocol "iGenomx_PCR". Afterward, secure DNA Plate back on the Magnetic Module. Resume protocol. (Press ENTER)')
            input("ARE YOU SURE? Press ENTER to continue.")
        else:
            self.protocol.pause('Resuspend each well in 50 ul Amplification "Reaction C" mix (primers added). Be sure to record which wells received which index primer. Execute thermocycler protocol "iGenomx_PCR". Afterward, secure DNA Plate back on the Magnetic Module. Resume protocol.')
        self.mag_mod.disengage()
        self.protocol.home()

        ### prepare for bead-based size selection
        self.set_flow_rate(48/4, 48, 48)
        self.mag_mod.engage()
        self.protocol.delay(minutes=8)
        self.p300.transfer(50, self.dna_plate.wells_by_name()['A2'].bottom(z=0.5), self.dna_plate.wells_by_name()['A3'])
        self.mag_mod.disengage()

    def Size_Selection(self):
        SSH_MODE = self.ssh_mode

        ### reagent locations
        Purification_Beads = self.post_pool_reagents.wells_by_name()['A8']
        Ethanol = [self.post_pool_reagents.wells_by_name()['A9'], self.post_pool_reagents.wells_by_name()['A10']]
        IDTE = self.post_pool_reagents.wells_by_name()['A12']

        ### first round of size selection (removes large reads)
        self.p300.pick_up_tip()
        self.set_flow_rate(48*3, 48*3, 48)
        self.p300.mix(10, 35, Purification_Beads)
        self.set_flow_rate(48, 48, 48)
        self.p300.aspirate(32.5, Purification_Beads)
        self.p300.dispense(32.5, self.dna_plate.wells_by_name()['A3'])
        self.set_flow_rate(48*3, 48*3, 48)
        self.p300.mix(10, 50, self.dna_plate.wells_by_name()['A3'])
        self.set_flow_rate(48, 48, 48)
        for i in range(5):
            self.p300.aspirate(50, self.dna_plate.wells_by_name()['A3'].bottom(z=1.0))
            self.p300.dispense(50, self.dna_plate.wells_by_name()['A3'].center())
        self.set_flow_rate(48*3, 48*3, 48)
        self.p300.mix(10, 50, self.dna_plate.wells_by_name()['A3'])
        self.p300.blow_out(self.dna_plate.wells_by_name()['A3'].top())
        self.p300.drop_tip()
        self.protocol.delay(minutes=10)
        self.mag_mod.engage(height=11.0)
        self.protocol.delay(minutes=6)
        self.set_flow_rate(48/4, 48/2, 48/2)
        self.p300.transfer(65, self.dna_plate.wells_by_name()['A3'], self.dna_plate.wells_by_name()['A4'], blow_out=True, blowout_location='destination well')
        self.set_flow_rate(1, 5, 5)
        self.p20.pick_up_tip()
        self.p20.aspirate(15, self.dna_plate.wells_by_name()['A3'].bottom())
        self.p20.dispense(10, self.dna_plate.wells_by_name()['A4'].top(z=-5.0))
        self.p20.dispense(5, self.dna_plate.wells_by_name()['A4'])
        self.p20.blow_out(self.dna_plate.wells_by_name()['A4'].top())
        self.p20.drop_tip()
        self.mag_mod.disengage()

        ### second round of size selection (removes small reads)
        self.p20.pick_up_tip()
        self.set_flow_rate(10, 40, 10)
        self.p20.mix(10, 10, Purification_Beads)
        self.set_flow_rate(10/2, 10/2, 10)
        self.p20.aspirate(10, Purification_Beads)
        self.p20.dispense(10, self.dna_plate.wells_by_name()['A4'])
        self.set_flow_rate(10*2, 10*2, 10)
        self.p20.mix(10, 20, self.dna_plate.wells_by_name()['A4'])
        self.set_flow_rate(10*2, 10*2, 10)
        for i in range(20):
            self.p20.aspirate(20, self.dna_plate.wells_by_name()['A4'].bottom(z=1.0))
            self.p20.dispense(20, self.dna_plate.wells_by_name()['A4'].top(z=-5.0))
        self.set_flow_rate(10*2, 10*2, 10)
        self.p20.mix(10, 20, self.dna_plate.wells_by_name()['A4'])
        self.p20.blow_out(self.dna_plate.wells_by_name()['A4'].top())
        self.p20.drop_tip()
        self.protocol.delay(minutes=10)
        self.mag_mod.engage()
        self.protocol.delay(minutes=6)
            # discard supernatant in empty column in case significant amount of beads are aspirated away
        self.set_flow_rate(48/4, 48/4, 48/4)
        self.p300.transfer(70, self.dna_plate.wells_by_name()['A4'], self.dna_plate.wells_by_name()['A7'], blow_out=True, blowout_location='destination well')

        ### wash with 80% ethanol twice
        self.set_flow_rate(48, 48, 48)
        self.p300.transfer(150, Ethanol[0], self.dna_plate.wells_by_name()['A4'])
        self.protocol.delay(seconds=30)
        self.discard_supernatant(150, self.dna_plate.wells_by_name()['A4'])
            # second time
        self.set_flow_rate(48, 48, 48)
        self.p300.transfer(150, Ethanol[1], self.dna_plate.wells_by_name()['A4'])
        self.protocol.delay(seconds=30)
        self.discard_supernatant(150, self.dna_plate.wells_by_name()['A4'])
            # guarantee discard of residual ethanol
        self.set_flow_rate(1, 48, 48)
        self.p20.pick_up_tip()
        self.p20.aspirate(15, self.dna_plate.wells_by_name()['A4'].bottom(z=-0.5))
        self.p20.drop_tip()
        self.set_flow_rate(48, 48, 48)
        self.protocol.delay(seconds=30)
        self.mag_mod.disengage()

        ### elute in IDTE
        self.set_flow_rate(48, 48, 48)
        self.p300.pick_up_tip()
        self.p300.aspirate(15, IDTE)
        self.p300.dispense(15, self.dna_plate.wells_by_name()['A4'])
        self.set_flow_rate(48*2, 48*6, 48)
        self.p300.mix(10, 20, self.dna_plate.wells_by_name()['A4'])
        self.p300.blow_out(self.dna_plate.wells_by_name()['A4'].top())
        self.p300.move_to(self.dna_plate.wells_by_name()['A4'].top(z=10.0))
        self.protocol.delay(minutes=5)
            # resuspend settled beads
        self.p300.mix(10, 20, self.dna_plate.wells_by_name()['A4'])
        self.p300.blow_out(self.dna_plate.wells_by_name()['A4'].top())
        self.p300.drop_tip()
        self.protocol.delay(minutes=5)
        if SSH_MODE:
            input('Resume to magstand final elution.')
            input("ARE YOU SURE? Press ENTER to continue.")
        else:
            self.protocol.pause('Resume to magstand final elution.')
        self.mag_mod.engage()
        self.protocol.delay(minutes=5)

        ### prepare for NGS prep
        self.set_flow_rate(12, 16, 16)
        self.p20.transfer(12, self.dna_plate.wells_by_name()['A4'].bottom(z=0.5), self.dna_plate.wells_by_name()['A12'])
        self.mag_mod.disengage()
        self.protocol.home()

    def set_flow_rate(self, aspirate=96, dispense=96, blow_out=96):
        self.p300.flow_rate.aspirate = aspirate
        self.p300.flow_rate.dispense = dispense
        self.p300.flow_rate.blow_out = blow_out
        self.p20.flow_rate.aspirate = aspirate
        self.p20.flow_rate.dispense = dispense
        self.p20.flow_rate.blow_out = blow_out

    def discard_supernatant(self, vol, well, trash=None):
        if trash is None:
            trash = self.protocol.fixed_trash['A1']
        if not self.p300.has_tip:
            self.p300.pick_up_tip()
        self.set_flow_rate(48/2, 48, 48)
        while vol > 0:
            if vol < 200:
                self.p300.aspirate(vol, well.bottom(z=0.5))
                self.p300.dispense(vol, trash)
            else:
                self.p300.aspirate(200, well.bottom(z=0.5))
                self.p300.dispense(200, trash)
            vol = vol - 200
        self.p300.drop_tip()
        self.set_flow_rate(48, 48, 48)

    def discard_residual(self, well):
        self.set_flow_rate(1, 48, 48)
        self.p20.pick_up_tip()
        self.p20.aspirate(5, well.bottom())
        self.set_flow_rate(1/2, 48, 48)
        self.p20.aspirate(2, well.bottom(z=-0.5))
        self.p20.drop_tip()
        self.set_flow_rate(48, 48, 48)
