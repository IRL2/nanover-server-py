import sys
from os import getcwd, path
from os.path import basename
import gradio as gr
from nanover.omni import OmniRunner
from nanover.omni.playback import PlaybackSimulation
from nanover.omni.openmm import OpenMMSimulation
from nanover.omni.record import record_from_server

imd_runner = None


def initialize_simulation_files(input_files, frame_interval, include_velocities, include_forces):
    simulation_files = []
    if input_files:
        input_files = input_files["files"]
        for file in input_files:
            if file.endswith(".xml"):
                simulation = OpenMMSimulation.from_xml_path(file, name=str(basename(file)))
                simulation.frame_interval = frame_interval
                simulation.include_velocities = include_velocities
                simulation.include_forces = include_forces
                simulation_files.append(simulation)
            else:
                raise ValueError(f"Invalid file type: {file.name}")
    return simulation_files


def initialize_recording_playbacks(trajectory_files, state_file):
    recording_playbacks = []
    trajectory_files = trajectory_files["files"] if trajectory_files else []
    state_file = state_file["files"] if state_file else []
    for i, traj_file in enumerate(trajectory_files):
        recording_playbacks.append(
            PlaybackSimulation(
                name=f"recording-playback_{i}",
                traj=traj_file,
                state=state_file[i] if i < len(state_file) else None,
            )
        )
    return recording_playbacks


def run_simulation(simulation_type, input_files, trajectory_files, state_file, server_name, port, simulation_fps,
                   frame_interval, start_paused, include_velocities, include_forces, record_trajectory,
                   trajectory_output_file, shared_state_file):
    global imd_runner
    try:
        simulation_files = initialize_simulation_files(input_files, frame_interval, include_velocities, include_forces)
        recording_playbacks = initialize_recording_playbacks(trajectory_files, state_file)

        imd_runner = OmniRunner.with_basic_server(
            *tuple(simulation_files + recording_playbacks), name=server_name, port=port
        )
        imd_runner.next()
        imd_runner.runner.play_step_interval = 1 / simulation_fps
        if start_paused:
            imd_runner.pause()

        if record_trajectory:
            if not trajectory_output_file or not shared_state_file:
                raise ValueError("Please provide both a trajectory output file and a shared state file.")
            record_from_server(
                f"localhost:{imd_runner.app_server.port}",
                trajectory_output_file,
                shared_state_file,
            )
        return f"Simulation started with type: settings: {locals()}"
    except Exception as e:
        if imd_runner:
            imd_runner.close()
        raise gr.Error(f"Error starting simulation: {e}")


def stop_simulation():
    global imd_runner
    try:
        if imd_runner:
            imd_runner.close()
        return "Simulation stopped!"
    except NameError:
        raise gr.Error("No simulation running!")


def toggle_visibility(choice):
    return gr.update(visible=choice == "From xml"), gr.update(visible=choice == "From recording")


def create_ui():
    with gr.Blocks(title='NanoVer', theme=gr.themes.Origin()) as demo:
        gr.Markdown("# NanoVer IMD Python Server GUI")

        simulation_type = gr.Radio(["From xml", "From recording"], label="Simulation Type", value="From xml")

        with gr.Row():
            with gr.Column(visible=True) as realtime_col:
                input_files = gr.MultimodalTextbox(label="Input Files (From xml)", file_count="multiple",
                                                   submit_btn=False)
            with gr.Column(visible=False) as playback_col:
                trajectory_files = gr.MultimodalTextbox(label="Trajectory Files (for playback)", file_count="multiple",
                                                        submit_btn=False)
                state_file = gr.MultimodalTextbox(label="State File (for playback)", file_count="multiple",
                                                  submit_btn=False)

        with gr.Row():
            with gr.Column():
                gr.Markdown("## Network")
                server_name = gr.Textbox(label="Server name", value="NanoVer-PY-GUI iMD Server")
                port = gr.Number(label="Port", value=38801)
                gr.Markdown("## Simulation")
                simulation_fps = gr.Slider(1, 60, value=30, label="Simulation FPS")
                frame_interval = gr.Number(label="Frame interval", value=5)
                start_paused = gr.Checkbox(label="Start simulation paused")
                include_velocities = gr.Checkbox(label="Include the velocities in the frames")
                include_forces = gr.Checkbox(label="Include the forces in the frames")

            with gr.Column():
                gr.Markdown("## Recording")
                with gr.Group():
                    record_trajectory = gr.Checkbox(label="Record Session")
                    trajectory_output_file = gr.Textbox(label="Trajectory output file path",
                                                        value=f"{path.join(getcwd(), 'filename.traj')}")
                    shared_state_file = gr.Textbox(label="Shared state file path",
                                                   value=f"{path.join(getcwd(), 'filename.state')}")

        run_button = gr.Button("Run the selected file!", variant='primary')
        stop_button = gr.Button("Stop the simulation!", variant='stop')
        output = gr.Textbox(label="Output")

        simulation_type.change(toggle_visibility, inputs=[simulation_type], outputs=[realtime_col, playback_col])
        run_button.click(run_simulation,
                         inputs=[simulation_type, input_files, trajectory_files, state_file, server_name, port,
                                 simulation_fps, frame_interval, start_paused, include_velocities, include_forces,
                                 record_trajectory, trajectory_output_file, shared_state_file], outputs=output)
        stop_button.click(stop_simulation, outputs=output)
    return demo


if __name__ == "__main__":
    demo = create_ui()
    demo.launch(inbrowser=True)