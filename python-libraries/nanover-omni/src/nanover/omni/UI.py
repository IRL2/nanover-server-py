import subprocess
import sys
# Function to check if a package is installed
def check_and_install(package):
    try:
        __import__(package)
    except ImportError:
        subprocess.check_call([sys.executable, "-m", "pip", "install", package])

# Check and install gradio
check_and_install("gradio")

import gradio as gr

def run_simulation(simulation_type, input_files, trajectory_files, state_file, verbosity, show_progression, server_name,
                   port, simulation_fps, frame_interval, force_interval, start_paused, include_velocities,
                   include_forces, record_stats, stats_file, stats_fps, record_trajectory,
                   trajectory_output_file, record_shared_state, shared_state_file):
    '''
    Runs a simulation with the given parameters.

    :param simulation_type: choose between "From xml" and "From recording"
    :param input_files: list of input files to run as live simulation
    :param trajectory_files: list of trajectory files to run as playback simulation
    :param state_file: list of state files to run as playback simulation
    :param verbosity: choose between "Normal", "Verbose", "Super verbose"
    :param show_progression: show simulation progression
    :param server_name: name of the server
    :param port: port number
    :param simulation_fps: frames per second
    :param frame_interval: frame interval
    :param force_interval: force interval
    :param start_paused: start simulation paused
    :param include_velocities: include velocities in the frames
    :param include_forces: include forces in the frames
    :param record_stats: record statistics
    :param stats_file: statistics file
    :param stats_fps: statistics frames per second
    :param record_trajectory: record trajectory
    :param trajectory_output_file: trajectory output file
    :param record_shared_state: record shared state
    :param shared_state_file: shared state file

    :return: a string with the simulation type and settings
    '''
    from nanover.omni import OmniRunner
    from nanover.omni.playback import PlaybackSimulation
    from nanover.omni.openmm import OpenMMSimulation
    from nanover.omni.record import record_from_server
    imd_runner = None
    global imd_runner

    # Initialize simulation files list
    simulation_files = list()
    for i in range(len(input_files)):
        if input_files[i].endswith('.xml'):
            # Create OpenMMSimulation from XML file
            simulation = OpenMMSimulation.from_xml_path(input_files[i], name=str(i))
            simulation.frame_interval = frame_interval
            simulation.include_velocities = include_velocities
            simulation.include_forces = include_forces
            simulation_files.append(simulation)
        else:
            return f"Invalid file type: {input_files[i].name}"

    # Initialize recording playbacks list
    recording_playbacks = list()
    trajectory_files = [] if trajectory_files is None else trajectory_files
    for i in range(len(trajectory_files)):
        # Create PlaybackSimulation from trajectory and state files
        recording_playbacks.append(PlaybackSimulation(name='recording-playback_'+str(i),
                                           traj=trajectory_files[i],
                                           state=state_file[i]))
    # Create OmniRunner with basic server
    imd_runner = OmniRunner.with_basic_server(*tuple(simulation_files+recording_playbacks), name=server_name, port=port)
    imd_runner.next()
    if start_paused:
        imd_runner.pause()
    if record_trajectory or record_shared_state:
        # Record from server if required
        record_from_server(f"localhost:{imd_runner.app_server.port}",trajectory_output_file, shared_state_file)

    return f"Simulation started with type: {simulation_type}, settings: {locals()}"

def stop_simulation():
    global imd_runner
    try:
        # Attempt to close the simulation runner
        imd_runner.close()
    except NameError as e:
        return e
    return "Simulation stopped!"
def create_ui():
    with gr.Blocks(theme=gr.themes.Soft(), fill_height=True) as demo:
        gr.Markdown("# Nanover IMD Python Server GUI")

        # Radio button to select simulation type
        simulation_type = gr.Radio(["From xml", "From recording"], label="Simulation Type", value="From xml")

        with gr.Row():
            with gr.Column(visible=True) as realtime_col:
                # File input for live simulation
                input_files = gr.File(label="Input Files (for From xml)", file_count="multiple")
            with gr.Column(visible=False) as playback_col:
                # File inputs for playback simulation
                trajectory_files = gr.File(label="Trajectory Files (for playback)", file_count="multiple")
                state_file = gr.File(label="State File (for playback)", file_count="multiple")

        with gr.Row():
            with gr.Column():
                gr.Markdown("## Verbosity")
                # Radio button for verbosity level
                verbosity = gr.Radio(["Normal", "Verbose", "Super verbose"], label="Verbosity Level", value="Normal")
                show_progression = gr.Checkbox(label="Show simulation progression")

                gr.Markdown("## Network")
                # Textbox for server name
                server_name = gr.Textbox(label="Server name", value="NanoVer-RS iMD Server")
                # Number input for port
                port = gr.Number(label="Port", value=38801)

                gr.Markdown("## Simulation")
                # Slider for simulation FPS
                simulation_fps = gr.Slider(1, 60, value=30, label="Simulation FPS")
                # Number inputs for frame and force intervals
                frame_interval = gr.Number(label="Frame interval", value=5)
                force_interval = gr.Number(label="Force interval", value=5)
                # Checkboxes for simulation options
                start_paused = gr.Checkbox(label="Start simulation paused")
                include_velocities = gr.Checkbox(label="Include the velocities in the frames")
                include_forces = gr.Checkbox(label="Include the forces in the frames")

            with gr.Column():
                gr.Markdown("## Recording")
                with gr.Group():
                    # Checkboxes and textboxes for recording options
                    record_stats = gr.Checkbox(label="Record statistics")
                    stats_file = gr.Textbox(label="Statistics file")
                    stats_fps = gr.Number(label="Statistics FPS", value=4)
                with gr.Group():
                    record_trajectory = gr.Checkbox(label="Record trajectory")
                    trajectory_output_file = gr.Textbox(label="Trajectory output file", value=f"{server_name.value}.traj")
                with gr.Group():
                    record_shared_state = gr.Checkbox(label="Record shared state")
                    shared_state_file = gr.Textbox(label="Shared state file", value=f"{server_name.value}.state")

        # Buttons to run and stop the simulation
        run_button = gr.Button("Run the selected file!")
        stop_button = gr.Button("Stop the simulation!")
        output = gr.Textbox(label="Output")

        def toggle_visibility(choice):
            # Toggle visibility of columns based on simulation type
            return gr.update(visible=choice == "From xml"), gr.update(visible=choice == "From recording")

        # Change event for simulation type radio button
        simulation_type.change(toggle_visibility, inputs=[simulation_type], outputs=[realtime_col, playback_col])

        # Click events for run and stop buttons
        run_button.click(
            run_simulation,
            inputs=[simulation_type, input_files, trajectory_files, state_file, verbosity, show_progression, server_name,
                    port, simulation_fps, frame_interval, force_interval, start_paused, include_velocities,
                    include_forces, record_stats, stats_file, stats_fps, record_trajectory,
                    trajectory_output_file, record_shared_state, shared_state_file],
            outputs=output
        )
        stop_button.click(stop_simulation, outputs=output)
    return demo
# Launch the Gradio interface
if __name__ == "__main__":
    ui = create_ui()
    ui.launch()