import streamlit as st
from streamlit_stl import stl_from_file
from utils.dataloader import log_event, log_submission
import time
def right_column():
    st.markdown('<h3 style="display: flex; align-items: center; justify-content: center; margin-top: -20px; text-align:center;">📸 Preview </h3>', unsafe_allow_html=True)

    if st.session_state['stl_path']:
        if "current_index" not in st.session_state:
            st.session_state.current_index = len(st.session_state['stl_files'])-1
        with col1:
            if st.session_state.current_index > 0:
                if st.button("←"):
                    st.session_state.current_index -= 1
                    st.session_state['switch'] = True

        with col3:
            if st.session_state.current_index < len(st.session_state['stl_files'])-1:
                if st.button("→"):
                    st.session_state.current_index += 1
                    st.session_state['switch'] = True
        with col2:
            path = st.session_state['stl_files'].keys()[st.session_state.current_index]
            if st.session_state['stl_files'] != {}:
                stl_from_file(path,st.session_state.get('stl_color', '#336fff'), 
                        auto_rotate=True, width=700, height=500,
                        cam_distance=100*(current_params['resolution']/50),cam_h_angle=45,cam_v_angle=75)
                schema = data["params_dict"].get(st.session_state['dict_key'], 1)
                st.session_state['slider_placeholder'].empty()
                with st.session_state['slider_placeholder']:
                    for param_key in schema:
                        labeled_slider(param_key,schema[param_key], st.session_state['current_params'],
                        preset=True)
         
            else:
              stl_from_file(st.session_state['stl_path'],st.session_state.get('stl_color', '#336fff'), 
                        auto_rotate=True, width=700, height=500,
                        cam_distance=100*(current_params['resolution']/50),cam_h_angle=45,cam_v_angle=75)
              
            
        # Download button
        try:
            design_ques_tab = st.columns([1, 1, 1])
            design_ques_list = ["design ques 1",]
            with design_ques_tab[1]:
                design_ques = st.selectbox("Design Question", design_ques_list, index=0)
                log_event(design_ques,'Chat mode')
            download_submit_tab = st.columns([1, 1, 1])
            with download_submit_tab[2]:
                if st.button('Submit'):
                    log_submission(st.session_state['confirmed_params'],design_ques, 'Pro Mode')
                    st.rerun() 
            with open(stl_path, 'rb') as f:
                
                with download_submit_tab[0]:
                    st.download_button('⬇️ Download STL', data=f.read(), file_name=stl_path, mime='model/stl')
                    log_event('Download','Pro mode')
        except Exception:
            st.warning('STL file not available for download. Generate again.')
    else:
        placeholder = st.empty()
        placeholder.markdown(
"""
<div style="
    width: 525px;
    height:350px;
    border: 2px dashed #ccc;
    display: flex;
    align-items: center;
    justify-content: center;
    margin: 0 auto;  /* center horizontally */
    text-align: center;
">
    Your STL will appear here
</div>
""",
unsafe_allow_html=True
)
