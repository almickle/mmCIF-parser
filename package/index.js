
export default function parse_mmCIF(text) {

    const ATOMS = text.split('\n').filter((line) => line.startsWith('ATOM')).map(array => array.split(' ').filter(element => element !== ''))
    const customLabels = ['', 'id', 'atom', 'atom_type', '', 'residue', 'chain', 'entity_index', 'residue_index', '', 'x', 'y', 'z', 'occupancy', 'isotropic_temperature_factor', 'formal_charge', 'author_residue_index', 'author_residue', 'author_chain', 'author_type', '' ]    
    const newAtoms = ATOMS.map((array) => {
        let object = {}
        array.forEach((entry, index) => {
            const label = customLabels[index]
            if(label !== '') {
                object = ({...object, [label]: entry})
            }
        })
        return (
            object
        )
    })

    return newAtoms
}
